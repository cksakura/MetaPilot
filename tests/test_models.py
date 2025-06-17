import pytest
import torch
import numpy as np
from src.models.pytorch_model import PeptideIntensityPredictor, PeptideDataset

def test_peptide_dataset():
    # Create dummy data
    features = np.random.rand(100, 15)
    targets = np.random.rand(100)
    weights = np.random.rand(100)
    
    # Create dataset
    dataset = PeptideDataset(features, targets, weights)
    
    # Test dataset length
    assert len(dataset) == 100
    
    # Test getitem
    sample = dataset[0]
    assert len(sample) == 3
    assert isinstance(sample[0], torch.Tensor)
    assert isinstance(sample[1], torch.Tensor)
    assert isinstance(sample[2], torch.Tensor)

def test_model_initialization():
    # Test model initialization
    model = PeptideIntensityPredictor(input_size=15)
    
    # Test model parameters
    assert isinstance(model, torch.nn.Module)
    assert hasattr(model, 'lstm1')
    assert hasattr(model, 'lstm2')
    assert hasattr(model, 'attention')
    assert hasattr(model, 'fc1')
    assert hasattr(model, 'fc2')
    assert hasattr(model, 'fc3')

def test_model_forward():
    # Create model
    model = PeptideIntensityPredictor(input_size=15)
    
    # Create dummy input
    x = torch.randn(32, 15)  # batch_size=32, features=15
    
    # Test forward pass
    output = model(x)
    assert output.shape == (32, 1)  # batch_size=32, output_size=1 