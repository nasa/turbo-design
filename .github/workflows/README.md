# GitHub Actions Workflows

## Documentation Build and Deploy

The `docs.yml` workflow automatically builds and deploys the Sphinx HTML documentation to GitHub Pages.

### Workflow Triggers

- **Push to main/master**: Builds and deploys documentation
- **Pull requests**: Builds documentation (no deployment) to verify changes
- **Manual dispatch**: Can be triggered manually from GitHub Actions tab

### Setup Requirements

To enable GitHub Pages deployment:

1. Go to your repository **Settings** → **Pages**
2. Under **Source**, select **GitHub Actions**
3. Save the settings

### Workflow Steps

1. **Build Job**
   - Checks out the repository
   - Sets up Python 3.12
   - Installs Sphinx and dependencies
   - Builds HTML documentation using `make html`
   - Uploads documentation as artifact

2. **Deploy Job** (only on push to main/master)
   - Deploys the built documentation to GitHub Pages
   - Documentation will be available at: `https://<username>.github.io/<repo-name>/`

### Local Documentation Build

To build documentation locally:

```bash
cd docs
make html
```

The built documentation will be in `docs/build/html/`. Open `docs/build/html/index.html` in a browser.

### Updating Documentation

Documentation is automatically regenerated when:
- Docstrings in Python files are updated
- RST files in `docs/` are modified
- `docs/conf.py` configuration is changed

The workflow uses mocked imports for heavy dependencies (Cantera, pyturbo) to speed up builds.
