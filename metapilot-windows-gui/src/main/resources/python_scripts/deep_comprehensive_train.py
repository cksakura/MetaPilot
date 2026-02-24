#!/usr/bin/env python3
"""
Deep Comprehensive two-head training script.

Reads the training seed TSV produced by DeepComprehensiveTrainingDatasetBuilder
(columns: peptide, length, miss_cleavage, gravy, basic_count, acidic_count, aromatic_count,
net_charge, detectability, n_genomes, n_proteins, n_species, n_genus, unique_species,
tax_sum_1..7, tax_max_1..7, tax_top2sum_1..7, fun_sum_1..5, fun_max_1..5, fun_top2sum_1..5,
fun_n_nodes_1..5, label_detect, label_intensity) and trains a two-head model:
- Detectability head: binary (observed vs not).
- Intensity head: regression (log1p(sum Normalised)).

Usage:
  python deep_comprehensive_train.py --input_tsv path/to/training_seed.tsv --output_dir path/to/training/model
  python deep_comprehensive_train.py --input_tsv training_seed.tsv --output_dir model --epochs 100 --batch_size 64

Requires: torch, pandas, numpy, scikit-learn
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

try:
    import torch
    import torch.nn as nn
    from torch.utils.data import DataLoader, TensorDataset
except ImportError:
    print("PyTorch is required. Install with: pip install torch", file=sys.stderr)
    sys.exit(1)

try:
    from sklearn.metrics import (
        average_precision_score,
        mean_absolute_error,
        r2_score,
        roc_auc_score,
    )
    from sklearn.model_selection import train_test_split
    from sklearn.preprocessing import StandardScaler
except ImportError:
    print("scikit-learn is required. Install with: pip install scikit-learn", file=sys.stderr)
    sys.exit(1)


# Feature columns in the same order as DeepComprehensiveTrainingDatasetBuilder TSV (excluding peptide and labels)
STATIC_COLUMNS = [
    "length", "miss_cleavage", "gravy", "basic_count", "acidic_count", "aromatic_count", "net_charge",
    "detectability",
]
CONTEXT_COLUMNS = [
    "n_genomes", "n_proteins", "n_species", "n_genus", "unique_species",
]
TAX_COLUMNS = []
for i in range(1, 8):
    TAX_COLUMNS.extend([f"tax_sum_{i}", f"tax_max_{i}", f"tax_top2sum_{i}"])
FUN_COLUMNS = []
for i in range(1, 6):
    FUN_COLUMNS.extend([f"fun_sum_{i}", f"fun_max_{i}", f"fun_top2sum_{i}", f"fun_n_nodes_{i}"])

FEATURE_COLUMNS = STATIC_COLUMNS + CONTEXT_COLUMNS + TAX_COLUMNS + FUN_COLUMNS
LABEL_DETECT = "label_detect"
LABEL_INTENSITY = "label_intensity"


def load_training_tsv(path: str | Path) -> pd.DataFrame:
    """Load training_seed.tsv and validate columns."""
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f"Training TSV not found: {path}")
    df = pd.read_csv(path, sep="\t", dtype={"peptide": str})
    missing = [c for c in FEATURE_COLUMNS + [LABEL_DETECT, LABEL_INTENSITY] if c not in df.columns]
    if missing:
        raise ValueError(f"Missing columns in {path}: {missing}")
    return df


def prepare_features_and_labels(
    df: pd.DataFrame,
    fill_detectability: str = "zero",
    val_frac: float = 0.15,
    test_frac: float = 0.1,
    random_state: int = 42,
):
    """
    Extract numeric features and labels, handle NaNs, split train/val/test.

    fill_detectability: 'zero' or 'median' for NaN detectability.
    Returns: X_train, X_val, X_test, y_detect_*, y_intensity_*, scaler, feature_names.
    """
    X = df[FEATURE_COLUMNS].copy()

    # Fill NaN detectability
    if X["detectability"].isna().any():
        if fill_detectability == "median":
            X["detectability"] = X["detectability"].fillna(X["detectability"].median())
        else:
            X["detectability"] = X["detectability"].fillna(0.0)

    # Inf/-inf can appear from log; clip or replace. Filling remaining NaNs with 0 is consistent
    # with the builder: for tax_*/fun_* columns, 0 means "no evidence" (absent taxa/functions).
    X = X.replace([np.inf, -np.inf], np.nan)
    if X.isna().any().any():
        X = X.fillna(0.0)

    y_detect = df[LABEL_DETECT].values.astype(np.float32)
    y_intensity = df[LABEL_INTENSITY].values.astype(np.float32)
    y_detect_bin = (y_detect > 0.5).astype(np.int32)

    # Train / val / test with stratified split on label_detect so positive ratio is similar in all splits
    idx = np.arange(len(X))
    try:
        idx_train, idx_rest = train_test_split(
            idx, test_size=val_frac + test_frac, random_state=random_state, stratify=y_detect_bin
        )
    except ValueError:
        idx_train, idx_rest = train_test_split(idx, test_size=val_frac + test_frac, random_state=random_state)
    val_ratio = test_frac / (val_frac + test_frac)
    try:
        idx_val, idx_test = train_test_split(
            idx_rest,
            test_size=val_ratio,
            random_state=random_state,
            stratify=y_detect_bin[idx_rest],
        )
    except ValueError:
        idx_val, idx_test = train_test_split(idx_rest, test_size=val_ratio, random_state=random_state)

    X_train, X_val, X_test = X.iloc[idx_train], X.iloc[idx_val], X.iloc[idx_test]
    y_d_train = y_detect[idx_train]
    y_d_val = y_detect[idx_val]
    y_d_test = y_detect[idx_test]
    y_i_train = y_intensity[idx_train]
    y_i_val = y_intensity[idx_val]
    y_i_test = y_intensity[idx_test]

    scaler = StandardScaler()
    X_train = scaler.fit_transform(X_train)
    X_val = scaler.transform(X_val)
    X_test = scaler.transform(X_test)

    return (
        X_train.astype(np.float32),
        X_val.astype(np.float32),
        X_test.astype(np.float32),
        y_d_train,
        y_d_val,
        y_d_test,
        y_i_train,
        y_i_val,
        y_i_test,
        scaler,
        FEATURE_COLUMNS,
    )


class TwoHeadModel(nn.Module):
    """
    Shared MLP backbone with two heads:
    - Detectability: binary classification (sigmoid).
    - Intensity: regression (identity).
    """

    def __init__(
        self,
        input_dim: int,
        hidden_dims: list[int] | None = None,
        dropout: float = 0.2,
    ):
        super().__init__()
        if hidden_dims is None:
            hidden_dims = [256, 128, 64]
        layers = []
        prev = input_dim
        for h in hidden_dims:
            layers.append(nn.Linear(prev, h))
            layers.append(nn.LayerNorm(h))
            layers.append(nn.ReLU(inplace=True))
            layers.append(nn.Dropout(dropout))
            prev = h
        self.backbone = nn.Sequential(*layers)
        self.head_detect = nn.Linear(prev, 1)
        self.head_intensity = nn.Linear(prev, 1)

    def forward(self, x: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
        z = self.backbone(x)
        logit_detect = self.head_detect(z)
        intensity = self.head_intensity(z).squeeze(-1)
        return logit_detect.squeeze(-1), intensity


def _masked_intensity_loss(
    pred_i: torch.Tensor,
    ib: torch.Tensor,
    pos_mask: torch.Tensor,
    mse: nn.MSELoss,
) -> torch.Tensor | None:
    """MSE on intensity only for positive (detected) samples. Returns None if no positives in batch."""
    if pos_mask.sum() == 0:
        return None
    return mse(pred_i[pos_mask], ib[pos_mask])


def train(
    X_train: np.ndarray,
    X_val: np.ndarray,
    y_d_train: np.ndarray,
    y_d_val: np.ndarray,
    y_i_train: np.ndarray,
    y_i_val: np.ndarray,
    device: torch.device,
    epochs: int = 100,
    batch_size: int = 64,
    lr: float = 1e-3,
    weight_detect: float = 1.0,
    weight_intensity: float = 1.0,
    early_stop_patience: int = 15,
) -> TwoHeadModel:
    """Train two-head model; returns model with best validation state loaded.
    Intensity loss is applied only to detected samples (label_detect==1); pseudo-negatives
    have no real intensity label and must not be used for the regression head.
    """
    input_dim = X_train.shape[1]
    model = TwoHeadModel(input_dim=input_dim).to(device)
    opt = torch.optim.AdamW(model.parameters(), lr=lr, weight_decay=1e-5)

    # pos_weight = neg_count / pos_count to handle class imbalance (negatives often >> positives)
    n_pos = int((y_d_train > 0.5).sum())
    n_neg = len(y_d_train) - n_pos
    pos_weight = torch.tensor([n_neg / max(1, n_pos)], dtype=torch.float32, device=device)
    bce = nn.BCEWithLogitsLoss(pos_weight=pos_weight)
    mse = nn.MSELoss()

    train_ds = TensorDataset(
        torch.from_numpy(X_train),
        torch.from_numpy(y_d_train),
        torch.from_numpy(y_i_train),
    )
    train_loader = DataLoader(train_ds, batch_size=batch_size, shuffle=True, drop_last=False)

    best_val_loss = float("inf")
    best_state: dict | None = None
    patience_counter = 0

    for epoch in range(epochs):
        model.train()
        train_loss = 0.0
        n_batches = 0
        for xb, db, ib in train_loader:
            xb, db, ib = xb.to(device), db.to(device), ib.to(device)
            pos_mask = db > 0.5

            opt.zero_grad()
            logit_d, pred_i = model(xb)
            loss_d = bce(logit_d, db)
            loss_i = _masked_intensity_loss(pred_i, ib, pos_mask, mse)
            if loss_i is not None:
                loss = weight_detect * loss_d + weight_intensity * loss_i
            else:
                loss = weight_detect * loss_d
            loss.backward()
            opt.step()
            train_loss += loss.item()
            n_batches += 1

        model.eval()
        with torch.no_grad():
            xv = torch.from_numpy(X_val).to(device)
            y_d_v = torch.from_numpy(y_d_val).to(device)
            y_i_v = torch.from_numpy(y_i_val).to(device)
            logit_d, pred_i = model(xv)
            val_loss_d = bce(logit_d, y_d_v).item()
            pos_val = y_d_v > 0.5
            if pos_val.sum() > 0:
                val_loss_i = mse(pred_i[pos_val], y_i_v[pos_val]).item()
                val_loss = weight_detect * val_loss_d + weight_intensity * val_loss_i
            else:
                val_loss = weight_detect * val_loss_d

        if (epoch + 1) % 10 == 0 or epoch == 0:
            avg_train = train_loss / max(1, n_batches)
            print(f"Epoch {epoch + 1}/{epochs}  train_loss={avg_train:.4f}  val_loss={val_loss:.4f}")

        if val_loss < best_val_loss:
            best_val_loss = val_loss
            best_state = {k: v.cpu().clone() for k, v in model.state_dict().items()}
            patience_counter = 0
        else:
            patience_counter += 1
            if early_stop_patience > 0 and patience_counter >= early_stop_patience:
                print(f"Early stopping at epoch {epoch + 1}")
                break

    if best_state is not None:
        model.load_state_dict(best_state)
    return model


def save_artifacts(
    output_dir: Path,
    model: TwoHeadModel,
    scaler: StandardScaler,
    feature_names: list[str],
    config: dict,
) -> None:
    """Save model state, scaler, feature list and config to output_dir."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    state_cpu = {k: v.cpu() for k, v in model.state_dict().items()}
    torch.save(state_cpu, output_dir / "model.pt")
    # Scaler: mean_ and scale_ are numpy arrays
    np.savez(
        output_dir / "scaler.npz",
        mean=scaler.mean_,
        scale=scaler.scale_,
        feature_names=np.array(feature_names, dtype=object),
    )
    with open(output_dir / "config.json", "w") as f:
        json.dump(config, f, indent=2)
    print(f"Saved model, scaler and config to {output_dir}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Train Deep Comprehensive two-head model (detectability + intensity).",
    )
    parser.add_argument("--input_tsv", required=True, help="Path to training_seed.tsv")
    parser.add_argument("--output_dir", required=True, help="Directory to save model.pt, scaler.npz, config.json")
    parser.add_argument("--epochs", type=int, default=100, help="Max training epochs")
    parser.add_argument("--batch_size", type=int, default=64, help="Batch size")
    parser.add_argument("--lr", type=float, default=1e-3, help="Learning rate")
    parser.add_argument("--val_frac", type=float, default=0.15, help="Validation fraction")
    parser.add_argument("--test_frac", type=float, default=0.1, help="Test fraction")
    parser.add_argument("--fill_detectability", choices=("zero", "median"), default="zero")
    parser.add_argument("--weight_detect", type=float, default=1.0, help="Loss weight for detectability head")
    parser.add_argument("--weight_intensity", type=float, default=1.0, help="Loss weight for intensity head")
    parser.add_argument("--early_stop_patience", type=int, default=15, help="Early stop after N epochs without improvement")
    parser.add_argument("--seed", type=int, default=42, help="Random seed")
    parser.add_argument("--device", default="cpu", help="Device: cpu or cuda")
    args = parser.parse_args()

    torch.manual_seed(args.seed)
    np.random.seed(args.seed)
    device = torch.device(args.device if torch.cuda.is_available() and "cuda" in args.device else "cpu")
    print(f"Using device: {device}")

    df = load_training_tsv(args.input_tsv)
    print(f"Loaded {len(df)} rows from {args.input_tsv}")

    (
        X_train, X_val, X_test,
        y_d_train, y_d_val, y_d_test,
        y_i_train, y_i_val, y_i_test,
        scaler,
        feature_names,
    ) = prepare_features_and_labels(
        df,
        fill_detectability=args.fill_detectability,
        val_frac=args.val_frac,
        test_frac=args.test_frac,
        random_state=args.seed,
    )
    print(f"Features: {len(feature_names)}  Train: {len(X_train)}  Val: {len(X_val)}  Test: {len(X_test)}")

    model = train(
        X_train, X_val,
        y_d_train, y_d_val,
        y_i_train, y_i_val,
        device=device,
        epochs=args.epochs,
        batch_size=args.batch_size,
        lr=args.lr,
        weight_detect=args.weight_detect,
        weight_intensity=args.weight_intensity,
        early_stop_patience=args.early_stop_patience,
    )

    config = {
        "input_tsv": str(Path(args.input_tsv).resolve()),
        "feature_columns": feature_names,
        "epochs": args.epochs,
        "batch_size": args.batch_size,
        "lr": args.lr,
        "seed": args.seed,
    }
    save_artifacts(
        Path(args.output_dir),
        model,
        scaler,
        feature_names,
        config,
    )

    # Test set metrics: AUC/PR-AUC for detectability (better than accuracy with imbalance); MSE/MAE/R² for intensity
    model.eval()
    with torch.no_grad():
        x_test = torch.from_numpy(X_test).to(device)
        logit_d, pred_i = model(x_test)
        prob_d = torch.sigmoid(logit_d).cpu().numpy()
        pred_i_np = pred_i.cpu().numpy()
    n_pos_test = int((y_d_test > 0.5).sum())
    print("--- Test set metrics ---")
    if n_pos_test > 0 and n_pos_test < len(y_d_test):
        try:
            auc = roc_auc_score(y_d_test, prob_d)
            pr_auc = average_precision_score(y_d_test, prob_d)
            print(f"  Detectability:  ROC-AUC={auc:.4f}  PR-AUC={pr_auc:.4f}")
        except Exception as e:
            print(f"  Detectability:  (AUC skipped: {e})")
    else:
        acc_d = (prob_d > 0.5).astype(np.float32) == y_d_test
        print(f"  Detectability:  accuracy={acc_d.mean():.4f}  (single-class set, AUC N/A)")
    mse_i = np.square(pred_i_np - y_i_test).mean()
    mae_i = mean_absolute_error(y_i_test, pred_i_np)
    r2_i = r2_score(y_i_test, pred_i_np)
    print(f"  Intensity (all): MSE={mse_i:.4f}  MAE={mae_i:.4f}  R²={r2_i:.4f}")
    if n_pos_test > 0:
        pos_test = y_d_test > 0.5
        mse_pos = np.square(pred_i_np[pos_test] - y_i_test[pos_test]).mean()
        mae_pos = mean_absolute_error(y_i_test[pos_test], pred_i_np[pos_test])
        r2_pos = r2_score(y_i_test[pos_test], pred_i_np[pos_test])
        print(f"  Intensity (pos): MSE={mse_pos:.4f}  MAE={mae_pos:.4f}  R²={r2_pos:.4f}")


if __name__ == "__main__":
    main()
