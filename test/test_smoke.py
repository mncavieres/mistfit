def test_import():
    import mistfit
    assert hasattr(mistfit, "fit_stars_with_minimint")

def test_cli_help():
    import subprocess, sys
    # Ensure the console script is installed in the CI venv
    result = subprocess.run([sys.executable, "-m", "mistfit.cli", "--help"], capture_output=True)
    assert result.returncode == 0
    assert b"Nested-sampling MIST + extinction" in result.stdout or b"usage:" in result.stdout
