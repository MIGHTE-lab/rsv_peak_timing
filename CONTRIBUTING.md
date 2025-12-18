# Contributing to RSV/Flu Timing Analysis

Thank you for your interest in contributing to this project! This document provides guidelines for contributing to the codebase.

## Code of Conduct

We are committed to providing a welcoming and inclusive environment. Please be respectful and professional in all interactions.

## How to Contribute

### Reporting Issues

If you find a bug or have a suggestion:

1. **Check existing issues** to avoid duplicates
2. **Open a new issue** with:
   - Clear, descriptive title
   - Detailed description of the problem/suggestion
   - Steps to reproduce (for bugs)
   - Expected vs actual behavior
   - Your environment (R version, OS, package versions)
   - Relevant code snippets or error messages

### Suggesting Enhancements

We welcome ideas for improvements:

1. **Open an issue** labeled "enhancement"
2. Describe the proposed feature clearly
3. Explain the use case and benefits
4. Consider implementation complexity

### Contributing Code

#### Before You Start

1. **Open an issue** to discuss major changes
2. **Fork the repository**
3. **Create a branch** for your work:
   ```bash
   git checkout -b feature/your-feature-name
   ```

#### Coding Standards

**R Code Style**:
- Follow [Tidyverse Style Guide](https://style.tidyverse.org/)
- Use meaningful variable names
- Add comments for complex logic
- Use `=` for assignment (not `<-`) for consistency with existing code

**Script Structure**:
```r
## -----------------------------------------------------------------------------
## Script name: script_name.R
##
## Purpose of script: Brief description
##
## Author: Your Name
##
## Date Created: YYYY-MM-DD
##
## Last Updated: YYYY-MM-DD
## Update notes: What changed
## -----------------------------------------------------------------------------

## 1. Load packages ----
library(tidyverse)

## 2. Set parameters ----
param1 = value1

## 3. Main analysis ----
# Clear section headers with ----

## 4. Save outputs ----
# ...
```

**Documentation**:
- Comment complex operations
- Use roxygen2-style documentation for functions:
  ```r
  #' Brief function description
  #'
  #' Detailed description if needed
  #'
  #' @param x Description of parameter x
  #' @param y Description of parameter y
  #' @return Description of return value
  #' @examples
  #' function_name(x = 1, y = 2)
  function_name = function(x, y) {
    # Function body
  }
  ```

**Testing**:
- Test your changes with multiple state-seasons
- Verify outputs match expected format
- Check edge cases (missing data, single observations, etc.)

#### Pull Request Process

1. **Update documentation** if needed:
   - README.md
   - CODE_GUIDE.md
   - Function documentation

2. **Test thoroughly**:
   - Run your modified scripts
   - Check for errors and warnings
   - Verify outputs are correct

3. **Commit changes**:
   ```bash
   git add file1.R file2.md
   git commit -m "Clear, descriptive commit message"
   ```

4. **Push to your fork**:
   ```bash
   git push origin feature/your-feature-name
   ```

5. **Open a Pull Request**:
   - Use a clear, descriptive title
   - Reference related issues (#123)
   - Describe changes in detail
   - List any breaking changes
   - Include example outputs if applicable

6. **Respond to feedback**:
   - Address reviewer comments
   - Make requested changes
   - Update PR description if scope changes

### Commit Message Guidelines

Use clear, imperative commit messages:

**Good**:
- "Add spatial correlation analysis"
- "Fix date parsing in EWS detection"
- "Update README with installation steps"
- "Refactor ridge regression function for clarity"

**Avoid**:
- "Fixed stuff"
- "Changes"
- "Update"

### Branch Naming

Use descriptive branch names:
- `feature/add-bootstrap-ci` - New features
- `fix/date-conversion-bug` - Bug fixes
- `docs/update-workflow` - Documentation
- `refactor/simplify-plotting` - Code refactoring

## Types of Contributions

### Documentation

Help improve documentation:
- Fix typos or unclear instructions
- Add examples
- Improve code comments
- Translate to other languages

### Code Quality

Improve code maintainability:
- Refactor for clarity
- Add error handling
- Optimize performance
- Add input validation

### New Features

Add new functionality:
- Additional statistical tests
- New visualization types
- Alternative modeling approaches
- Enhanced data processing

### Bug Fixes

Fix identified issues:
- Correct calculation errors
- Handle edge cases
- Fix plotting issues
- Resolve package conflicts

## Development Setup

1. **Clone your fork**:
   ```bash
   git clone https://github.com/your-username/rsv_peak_timing.git
   cd rsv_peak_timing
   ```

2. **Add upstream remote**:
   ```bash
   git remote add upstream https://github.com/MIGHTE-lab/rsv_peak_timing.git
   ```

3. **Keep your fork updated**:
   ```bash
   git fetch upstream
   git merge upstream/main
   ```

4. **Install dependencies**:
   ```r
   source("docs/INSTALLATION.md")  # Follow installation guide
   ```

5. **Run tests** (if available):
   ```r
   source("R/utils/test_installation.R")
   ```

## Code Review Process

All contributions are reviewed before merging:

1. **Automated checks** (if set up):
   - Code style
   - Basic tests

2. **Maintainer review**:
   - Code quality
   - Documentation
   - Compatibility
   - Scientific validity

3. **Feedback and revisions**:
   - Address comments
   - Make requested changes
   - Discuss disagreements respectfully

4. **Approval and merge**:
   - Once approved, maintainers will merge
   - Your contribution will be acknowledged

## Recognition

Contributors will be acknowledged in:
- CITATION.cff file
- Repository contributors list
- Publication acknowledgments (for significant contributions)

## Questions?

- **General questions**: Open an issue with the "question" label
- **Direct contact**: g.dewey@northeastern.edu
- **Discussions**: Use GitHub Discussions (if enabled)

## License

By contributing, you agree that your contributions will be licensed under the MIT License.

---

Thank you for contributing to this project! Your efforts help improve public health surveillance and epidemic preparedness. 🎉
