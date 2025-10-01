# Contributing to SparseQEEcNEO.jl

Thank you for your interest in contributing to SparseQEEcNEO.jl! This document provides guidelines for contributing to the project.

## Table of Contents

- [Code of Conduct](#code-of-conduct)
- [Getting Started](#getting-started)
- [Development Workflow](#development-workflow)
- [Code Style Guidelines](#code-style-guidelines)
- [Testing Requirements](#testing-requirements)
- [Documentation Standards](#documentation-standards)
- [Submitting Changes](#submitting-changes)

## Code of Conduct

We are committed to providing a welcoming and inclusive environment for all contributors. Please be respectful and constructive in all interactions.

## Getting Started

### Prerequisites

1. **Julia 1.6+**: Install from [julialang.org](https://julialang.org/)
2. **PySCF with NEO support**: Install from [https://github.com/theorychemyang/pyscf](https://github.com/theorychemyang/pyscf)
3. **Git**: For version control

### Environment Setup

```bash
# Clone the repository
git clone https://github.com/fzahari/SparseQEEcNEO.git
cd SparseQEEcNEO

# Install PySCF with NEO support
git clone https://github.com/theorychemyang/pyscf.git
cd pyscf && pip install -e . && cd ..

# Set up Julia environment
julia --project=.
julia> using Pkg; Pkg.instantiate()

# Set environment variables
source setup_env.sh

# Verify installation
julia scripts/check_dependencies.jl
```

## Development Workflow

1. **Fork** the repository on GitHub
2. **Create** a feature branch: `git checkout -b feature/your-feature-name`
3. **Develop** your changes following our code style guidelines
4. **Test** your changes thoroughly
5. **Document** your changes
6. **Submit** a pull request

### Branch Naming Conventions

- `feature/description` - New features
- `bugfix/description` - Bug fixes
- `docs/description` - Documentation improvements
- `refactor/description` - Code refactoring
- `test/description` - Test improvements

## Code Style Guidelines

We follow **Clean Code principles** throughout the codebase:

### Function Guidelines

- **Single Responsibility**: Each function should do one thing well
- **Small Functions**: Keep functions ≤ 20 lines when possible
- **Descriptive Names**: Use clear, self-documenting function names
- **Few Parameters**: Aim for ≤ 3 parameters; use structs for complex data

### Naming Conventions

```julia
# Good
function calculate_nuclear_position_expectation_value(nuclear_mean_field, origin)
const DEFAULT_CONVERGENCE_THRESHOLD = 1e-6
struct CNEOCalculation

# Avoid
function calc_pos(mf, o)  # Too abbreviated
const THR = 1e-6  # Magic number without context
struct Calc  # Too generic
```

### Documentation Requirements

Every function must have a docstring:

```julia
"""
    function_name(param1, param2)

Brief description of what the function does.

# Arguments
- `param1::Type`: Description of parameter 1
- `param2::Type`: Description of parameter 2

# Returns
- `ReturnType`: Description of what is returned

# Examples
```julia
result = function_name("example", 42)
```
"""
function function_name(param1, param2)
    # Implementation
end
```

### Constants and Configuration

- Use named constants instead of magic numbers
- Group related constants together
- Use ALL_CAPS for constants

```julia
# Clean Code Constants
const DEFAULT_CNEO_CONVERGENCE_THRESHOLD = 1e-6
const DEFAULT_CNEO_MAX_ITERATIONS = 50
const POSITION_DIMENSION = 3
```

## Testing Requirements

### Test Coverage

- **All new functions** must have corresponding tests
- **Integration tests** for workflows involving multiple components
- **Edge case testing** for error conditions and boundary values

### Test Organization

```
test/
├── runtests.jl                    # Main test runner
├── test_[module_name].jl          # Unit tests for each module
└── integration/                   # Integration tests
    ├── test_cneo_qee_workflow.jl
    └── test_quantum_export.jl
```

### Running Tests

```bash
# Run all tests
julia test/runtests.jl

# Run specific test module
julia test/test_cneo_methods.jl

# Run with coverage
julia --code-coverage=user test/runtests.jl
```

## Documentation Standards

### API Documentation

- Use clear, concise language
- Include examples for complex functions
- Document all parameters and return values
- Explain the scientific/mathematical background when relevant

### User Guides

- Provide step-by-step tutorials
- Include complete, runnable examples
- Explain the scientific context
- Show expected outputs

### Code Comments

- Explain **why**, not **what**
- Use comments sparingly for self-documenting code
- Add comments for complex algorithms or scientific formulations

## Submitting Changes

### Pull Request Guidelines

1. **Description**: Provide a clear description of changes
2. **Tests**: Include tests for new functionality
3. **Documentation**: Update relevant documentation
4. **Clean Code**: Ensure code follows our style guidelines
5. **No Breaking Changes**: Maintain backward compatibility when possible

### Pull Request Template

```markdown
## Description
Brief description of changes and motivation.

## Type of Change
- [ ] Bug fix
- [ ] New feature
- [ ] Documentation update
- [ ] Performance improvement
- [ ] Code refactoring

## Testing
- [ ] All existing tests pass
- [ ] New tests added for new functionality
- [ ] Integration tests updated if needed

## Documentation
- [ ] Code is self-documenting with clear function names
- [ ] Docstrings added/updated for public functions
- [ ] User documentation updated if needed

## Clean Code Compliance
- [ ] Functions are small and focused (≤ 20 lines preferred)
- [ ] Descriptive naming throughout
- [ ] Magic numbers replaced with named constants
- [ ] Single responsibility principle followed

## Checklist
- [ ] Code follows project style guidelines
- [ ] Self-review completed
- [ ] Comments explain complex logic
- [ ] No duplicate code
```

### Review Process

All pull requests require:
- **Code review** by maintainers
- **Passing tests** on CI
- **Clean Code compliance** verification
- **Documentation** completeness check

## Development Tools

### Quality Assurance

```bash
# Check dependencies
julia scripts/check_dependencies.jl

# Run benchmarks
julia scripts/benchmark_suite.jl

# Analyze code quality
julia scripts/analyze_clean_code.jl
```

### Debugging

- Use descriptive error messages
- Include context in error handling
- Test error conditions explicitly

## Questions and Support

- **Issues**: Use GitHub Issues for bugs and feature requests
- **Discussions**: Use GitHub Discussions for questions
- **Documentation**: Check the README.md and WARP.md files first

## Recognition

Contributors are recognized in:
- `CITATION.bib` file
- Project documentation
- Release notes

Thank you for contributing to SparseQEEcNEO.jl! Your contributions help advance quantum chemistry and quantum computing research.