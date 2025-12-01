# DNAbyte Testing Guide

## Running Tests

### Install test dependencies

```bash
# With uv
uv sync --extra test

# Or install dev dependencies
uv sync --extra dev
```

### Run all tests

```bash
# With uv
uv run pytest

# Or if pytest is installed
pytest
```

### Run specific test categories

```bash
# Run only unit tests
uv run pytest tests/unit/

# Run only integration tests
uv run pytest tests/integration/

# Run only fast tests (exclude slow tests)
uv run pytest -m "not slow"

# Run tests in parallel (faster)
uv run pytest -n auto
```

### Run with coverage

```bash
# Generate coverage report
uv run pytest --cov=dnabyte --cov-report=html

# View coverage report
open htmlcov/index.html
```

### Run specific test file or function

```bash
# Run specific test file
uv run pytest tests/unit/test_binarize.py

# Run specific test function
uv run pytest tests/unit/test_binarize.py::TestBinarize::test_binarize_text

# Run tests matching a pattern
uv run pytest -k "binarize"
```

## Test Structure

```
tests/
├── conftest.py              # Shared fixtures and configuration
├── unit/                    # Fast, isolated unit tests
│   ├── test_binarize.py
│   ├── test_encode.py
│   └── test_data_classes.py
├── integration/             # Tests combining multiple components
│   └── test_encoding_pipeline.py
```

## Writing Tests

### Example Unit Test

```python
import pytest
from dnabyte.data_classes.base import Data

def test_data_creation():
    """Test creating a Data object."""
    data = Data("test", data_type="text")
    assert data is not None
    assert data.data_type == "text"
```

### Example Integration Test

```python
@pytest.mark.integration
def test_encode_decode_roundtrip(sample_text_data):
    """Test that encode -> decode returns original data."""
    # Binarize
    data = Data(sample_text_data, data_type="text")
    binarizer = Binarize()
    binary = binarizer.binarize(data)
    
    # Encode
    encoder = Encode(method="linear_chain")
    encoded = encoder.encode(binary)
    
    # Decode
    decoded = encoder.decode(encoded)
    
    # Verify
    assert decoded == binary
```

### Example Parametrized Test

```python
@pytest.mark.parametrize("encoding_method", [
    "linear_chain",
    "poly_chain",
    "max_density",
])
def test_different_encodings(encoding_method, sample_data):
    """Test different encoding methods."""
    encoder = Encode(method=encoding_method)
    result = encoder.encode(sample_data)
    assert result is not None
```

## Best Practices

1. **Test at multiple levels**: Unit, integration, and end-to-end
2. **Use fixtures**: Share common test data and setup
3. **Mark slow tests**: Use `@pytest.mark.slow` for long-running tests
4. **Parametrize**: Test multiple scenarios with one test function
5. **Assert roundtrips**: Encode -> decode should return original data
6. **Test edge cases**: Empty data, large data, special characters
7. **Keep tests fast**: Unit tests should run in milliseconds
8. **Use meaningful names**: Test names should describe what they test

## Continuous Integration

Tests run automatically on GitHub Actions for every push and pull request.
See `.github/workflows/test.yml` for configuration.
