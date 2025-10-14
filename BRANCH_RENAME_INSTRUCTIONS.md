# Branch Renaming Instructions

This repository is in the process of renaming the `master` branch to `main` to use more inclusive terminology.

## Steps to Complete the Branch Rename

These steps must be performed by a repository administrator with appropriate GitHub permissions:

### 1. Rename the master branch to main on GitHub

**Option A: Using GitHub Web Interface**
1. Go to the repository on GitHub: https://github.com/cdshaffer/starterator
2. Navigate to Settings → Branches
3. Click the pencil icon next to the `master` branch
4. Rename it to `main`
5. Click "Rename branch"
6. GitHub will automatically update pull requests and branch protections

**Option B: Using Git Command Line**
```bash
# Clone the repository if you haven't already
git clone https://github.com/cdshaffer/starterator.git
cd starterator

# Rename the local master branch to main
git branch -m master main

# Push the new main branch to remote
git push -u origin main

# Set main as the default branch on GitHub (via Settings → Branches)
# Then delete the old master branch
git push origin --delete master
```

### 2. Update the default branch on GitHub
1. Go to Settings → Branches
2. Change the default branch from `master` to `main`
3. Confirm the change

### 3. Update local clones (for all contributors)
Anyone who has cloned this repository should update their local copy:

```bash
# Fetch the latest changes
git fetch origin

# Switch to the new main branch
git checkout main

# Set up tracking
git branch -u origin/main main

# Optionally, delete the old local master branch
git branch -d master
```

### 4. Update any CI/CD pipelines or scripts
Review and update any continuous integration configurations, deployment scripts, or other automation that references the `master` branch.

## What This PR Does

This pull request prepares the repository for the branch rename by:
- Providing these instructions for completing the rename
- Ensuring no hardcoded references to `master` exist in the codebase (already verified - none found)

## Additional Notes

- The repository analysis shows no hardcoded references to `master` branch in code files
- All git operations should work seamlessly after the rename
- GitHub will automatically redirect most references from `master` to `main`
