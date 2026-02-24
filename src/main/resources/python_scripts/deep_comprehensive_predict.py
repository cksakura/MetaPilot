#!/usr/bin/env python3
"""
Run the Deep Comprehensive two-head model on a predict_dataset.tsv and write predictions.

Reads the same feature columns as training (no labels). Outputs TSV: peptide, prob_detect, pred_intensity.

Usage:
  python deep_comprehensive_predict.py --model_dir path/to/training/model --predict_tsv path/to/predict_dataset.tsv --output_tsv path/to/predictions.tsv

Requires: torch, pandas, numpy (same as deep_comprehensive_train).
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

try:
    import torch
    import torch.nn as nn
except ImportError:
    print("PyTorch is required. Install with: pip install torch", file=sys.stderr)
    sys.exit(1)

# Reuse feature columns and model definition from training script
try:
    from deep_comprehensive_train import FEATURE_COLUMNS, TwoHeadModel
except ImportError:
    # Fallback if not run from same directory
    FEATURE_COLUMNS = None
    TwoHeadModel = None


def _ensure_feature_columns():
    if FEATURE_COLUMNS is not None:
        return FEATURE_COLUMNS
    # Inline same list as in deep_comprehensive_train
    static = [
        "length", "miss_cleavage", "gravy", "basic_count", "acidic_count", "aromatic_count", "net_charge",
        "detectability",
    ]
    context = ["n_genomes", "n_proteins", "n_species", "n_genus", "unique_species"]
    tax = []
    for i in range(1, 8):
        tax.extend([f"tax_sum_{i}", f"tax_max_{i}", f"tax_top2sum_{i}"])
    fun = []
    for i in range(1, 6):
        fun.extend([f"fun_sum_{i}", f"fun_max_{i}", f"fun_top2sum_{i}", f"fun_n_nodes_{i}"])
    return static + context + tax + fun


def _ensure_model_class():
    # Use globals() so we don't shadow with the inline class name (avoid UnboundLocalError)
    model_cls = globals().get("TwoHeadModel")
    if model_cls is not None:
        return model_cls
    import torch.nn as nn

    class _TwoHeadModelFallback(nn.Module):
        def __init__(self, input_dim, hidden_dims=None, dropout=0.2):
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

        def forward(self, x):
            z = self.backbone(x)
            logit_d = self.head_detect(z).squeeze(-1)
            intensity = self.head_intensity(z).squeeze(-1)
            return logit_d, intensity

    return _TwoHeadModelFallback


def load_scaler(model_dir: Path):
    npz = np.load(model_dir / "scaler.npz", allow_pickle=True)
    mean = np.asarray(npz["mean"], dtype=np.float32)
    scale = np.asarray(npz["scale"], dtype=np.float32)
    return mean, scale


def main() -> None:
    parser = argparse.ArgumentParser(description="Run Deep Comprehensive model on predict_dataset.tsv")
    parser.add_argument("--model_dir", required=True, help="Directory with model.pt and scaler.npz")
    parser.add_argument("--predict_tsv", required=True, help="Path to predict_dataset.tsv")
    parser.add_argument("--output_tsv", required=True, help="Path to write predictions.tsv")
    parser.add_argument("--device", default="cpu", help="Device: cpu or cuda")
    args = parser.parse_args()

    model_dir = Path(args.model_dir)
    if not model_dir.is_dir():
        print(f"Model dir not found: {model_dir}", file=sys.stderr)
        sys.exit(1)
    model_path = model_dir / "model.pt"
    if not model_path.is_file():
        print(f"model.pt not found in {model_dir}", file=sys.stderr)
        sys.exit(1)
    scaler_path = model_dir / "scaler.npz"
    if not scaler_path.is_file():
        print(f"scaler.npz not found in {model_dir}", file=sys.stderr)
        sys.exit(1)

    predict_tsv = Path(args.predict_tsv)
    if not predict_tsv.is_file():
        print(f"Predict TSV not found: {predict_tsv}", file=sys.stderr)
        sys.exit(1)

    feature_columns = _ensure_feature_columns()
    TwoHeadModelCls = _ensure_model_class()

    # Use CPU if user asked for cpu, or if CUDA is not available / not compatible (e.g. sm_120 unsupported)
    use_cuda = args.device.strip().lower() == "cuda" and torch.cuda.is_available()
    if use_cuda:
        try:
            torch.zeros(1).cuda()
        except Exception:
            use_cuda = False
    device = torch.device("cuda" if use_cuda else "cpu")
    mean, scale = load_scaler(model_dir)
    df = pd.read_csv(predict_tsv, sep="\t", dtype={"peptide": str})
    missing = [c for c in feature_columns if c not in df.columns]
    if missing:
        print(f"Missing columns in predict TSV: {missing}", file=sys.stderr)
        sys.exit(1)
    X = df[feature_columns].fillna(0.0).replace([np.inf, -np.inf], np.nan).fillna(0.0).astype(np.float32)
    X = ((X - mean) / scale).astype(np.float32)
    peptides = df["peptide"].values

    input_dim = len(feature_columns)
    model = TwoHeadModelCls(input_dim=input_dim)
    state = torch.load(model_path, map_location="cpu")
    model.load_state_dict(state)
    model.to(device)
    model.eval()

    with torch.no_grad():
        x = torch.from_numpy(X.values).to(device=device, dtype=torch.float32)
        logit_d, pred_i = model(x)
        prob_d = torch.sigmoid(logit_d).cpu().numpy()
        pred_i_np = pred_i.cpu().numpy()

    out = Path(args.output_tsv)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w") as f:
        f.write("peptide\tprob_detect\tpred_intensity\n")
        for i in range(len(peptides)):
            f.write(f"{peptides[i]}\t{prob_d[i]:.6f}\t{pred_i_np[i]:.6f}\n")
    print(f"Wrote {len(peptides)} predictions to {out}")


if __name__ == "__main__":
    main()
