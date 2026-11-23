# Contributing to this project

Thank you for your interest in contributing! This document explains how to propose changes.

## Fork-based workflow (required)

**All contributions must be made from a fork of this repository.** Direct pushes to this repository, and branches created on it by contributors, are not accepted. Only maintainers have write access to the main repository; everyone else (including regular contributors) works through a fork and a pull request.

### 1. Fork the repository

Click **Fork** at the top right of the repository page on GitHub. This creates a copy under your own account.

### 2. Clone your fork

```bash
git clone https://github.com/<your-username>/<repository>.git
cd <repository>
```

### 3. Add the original repository as `upstream`

```bash
git remote add upstream https://github.com/<owner>/<repository>.git
git fetch upstream
```

### 4. Create a branch for your change

Never work directly on `main`. Create a dedicated branch from an up-to-date `main`:

```bash
git checkout main
git pull upstream main
git checkout -b <type>/<short-description>
```

Suggested branch prefixes: `feature/`, `fix/`, `docs/`, `refactor/`, `test/`.

### 5. Make your changes

- Keep each pull request focused on a single topic.
- Follow the existing code style and conventions of the project.
- Add or update tests and documentation when relevant.
- Make sure the build and the test suite pass locally before submitting.

### 6. Commit

Write clear, descriptive commit messages (imperative mood, a short summary line of at most ~72 characters, followed by a body if needed):

```bash
git add <files>
git commit -m "Short summary of the change"
```

### 7. Keep your fork in sync

Before opening your pull request, rebase on the latest upstream `main`:

```bash
git fetch upstream
git rebase upstream/main
```

### 8. Push to your fork

```bash
git push origin <type>/<short-description>
```

### 9. Open a pull request

On GitHub, open a pull request from your fork's branch to the `main` branch of this repository. In the description:

- Explain **what** the change does and **why** it is needed.
- Link related issues (e.g. `Closes #123`).
- Mention any breaking change, new dependency or limitation.

## Review process

- A maintainer will review your pull request; please allow a few days for a first response.
- Address review comments by pushing additional commits to the same branch of your fork. The pull request updates automatically.
- Once approved, a maintainer will merge it. Maintainers may ask you to squash or rebase commits beforehand.

## Reporting bugs and requesting features

Before starting significant work, please **open an issue** to discuss it, so that effort is not spent on a change that cannot be accepted. When reporting a bug, include:

- the version or commit used,
- your environment (OS, compiler/interpreter, dependency versions),
- minimal steps to reproduce,
- expected versus observed behaviour.


## License

By contributing, you agree that your contributions will be licensed under the same license as this project (see `LICENSE`).
