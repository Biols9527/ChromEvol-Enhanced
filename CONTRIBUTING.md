# 🤝 Contributing Guidelines

We welcome contributions to the ChromEvol-Enhanced Chromosome Evolution Analysis project! This document provides guidelines for contributing.

## 🚀 Getting Started

### Prerequisites
- Python 3.7+
- Git
- Basic understanding of phylogenetics and chromosome evolution

### Setting Up Development Environment
```bash
# Fork and clone the repository
git clone https://github.com/your-username/chromevol-enhanced-analysis.git
cd chromevol-enhanced-analysis

# Create a virtual environment
python -m venv dev_env
source dev_env/bin/activate  # On Windows: dev_env\Scripts\activate

# Install development dependencies
pip install -r requirements.txt
pip install -r requirements-dev.txt  # Additional dev tools
```

## 📋 Types of Contributions

### 🐛 Bug Reports
- Use the GitHub issue tracker
- Provide detailed description and reproduction steps
- Include system information and error messages

### 💡 Feature Requests
- Discuss new features in GitHub issues first
- Provide clear use cases and requirements
- Consider implementation feasibility

### 🔧 Code Contributions
- Bug fixes
- New features
- Performance improvements
- Documentation updates

### 📚 Documentation
- API documentation
- User guides
- Tutorials and examples
- Code comments

## 🔄 Development Workflow

### 1. Create a Branch
```bash
git checkout -b feature/your-feature-name
# or
git checkout -b fix/bug-description
```

### 2. Make Changes
- Follow code style guidelines
- Add tests for new functionality
- Update documentation as needed

### 3. Test Your Changes
```bash
# Run basic tests
python src/ancestral_reconstruction.py --help

# Test with example data
python src/ancestral_reconstruction.py \
    --tree data/pruned_phylogeny.nwk \
    --counts data/chromosome_counts.csv \
    --out_image test_output.svg
```

### 4. Commit Changes
```bash
git add .
git commit -m "feat: add new chromosome visualization feature"
```

Use conventional commit messages:
- `feat:` new features
- `fix:` bug fixes
- `docs:` documentation changes
- `style:` formatting changes
- `refactor:` code refactoring
- `test:` adding tests
- `chore:` maintenance tasks

### 5. Push and Create Pull Request
```bash
git push origin feature/your-feature-name
```

Then create a pull request on GitHub.

## 📝 Code Style Guidelines

### Python Code Style
- Follow PEP 8
- Use meaningful variable names
- Add docstrings for functions and classes
- Maximum line length: 88 characters

### Example:
```python
def calculate_transition_probability(rate_matrix: np.ndarray, 
                                   branch_length: float) -> np.ndarray:
    """
    Calculate transition probabilities using matrix exponentiation.
    
    Args:
        rate_matrix: Q-matrix of evolutionary rates
        branch_length: Phylogenetic branch length
        
    Returns:
        Transition probability matrix P(t) = exp(Q*t)
    """
    return scipy.linalg.expm(rate_matrix * branch_length)
```

### Documentation Style
- Use clear, concise language
- Include examples for complex concepts
- Update documentation when changing functionality

## 🧪 Testing Guidelines

### Manual Testing
- Test with different input data formats
- Verify output file generation
- Check visualization rendering

### Test Cases to Consider
- Edge cases (very small/large chromosome numbers)
- Missing data handling
- Invalid input formats
- Performance with large datasets

## 📊 Performance Considerations

### Optimization Principles
- Use vectorized operations (NumPy/SciPy)
- Minimize memory allocation in loops
- Cache expensive computations
- Profile code for bottlenecks

### Memory Management
- Handle large datasets efficiently
- Use generators for data streaming
- Clear unused variables in long-running functions

## 🔬 Scientific Accuracy

### Biological Validity
- Ensure methods are scientifically sound
- Validate against known test cases
- Consider edge cases in biological data

### Statistical Rigor
- Implement proper statistical tests
- Handle uncertainty appropriately
- Document assumptions and limitations

## 📚 Documentation Standards

### Code Documentation
```python
class ChromEvolutionModel:
    """
    ChromEvol-inspired chromosome evolution model.
    
    This class implements maximum likelihood estimation of chromosome
    evolution parameters including gain, loss, and duplication rates.
    
    Attributes:
        max_chr: Maximum chromosome number to consider
        params: Dictionary of evolutionary rate parameters
        
    Example:
        >>> model = ChromEvolutionModel(max_chromosome_number=50)
        >>> model.set_parameters(gain=0.1, loss=0.2)
    """
```

### User Documentation
- Clear step-by-step instructions
- Practical examples
- Troubleshooting guides
- FAQ sections

## 🐛 Reporting Issues

### Bug Report Template
```markdown
**Bug Description**
A clear description of the bug.

**To Reproduce**
Steps to reproduce the behavior:
1. Run command '...'
2. With input file '...'
3. Error occurs

**Expected Behavior**
What you expected to happen.

**Environment**
- OS: [e.g., macOS 12.0]
- Python version: [e.g., 3.9.7]
- Package versions: [run `pip freeze`]

**Additional Context**
Any other context about the problem.
```

### Feature Request Template
```markdown
**Feature Description**
A clear description of the desired feature.

**Use Case**
Explain why this feature would be useful.

**Proposed Implementation**
If you have ideas about implementation.

**Alternatives Considered**
Other approaches you've considered.
```

## 🔄 Review Process

### Pull Request Checklist
- [ ] Code follows style guidelines
- [ ] Tests pass
- [ ] Documentation updated
- [ ] Commit messages are clear
- [ ] No merge conflicts

### Review Criteria
- Code quality and readability
- Scientific accuracy
- Performance impact
- Documentation completeness
- Test coverage

## 🏆 Recognition

Contributors will be:
- Listed in the CONTRIBUTORS.md file
- Mentioned in release notes
- Credited in academic papers (for significant contributions)

## 📞 Getting Help

### Communication Channels
- **GitHub Issues**: Bug reports and feature requests
- **GitHub Discussions**: General questions and ideas
- **Email**: [maintainer-email] for sensitive issues

### Development Questions
- Check existing issues and documentation first
- Provide context and specific questions
- Include relevant code snippets

## 📜 License

By contributing, you agree that your contributions will be licensed under the MIT License.

---

Thank you for contributing to the ChromEvol-Enhanced Analysis project! 🧬✨
