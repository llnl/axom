import pycore


def test_about_output(capfd):
    pycore.about()
    captured = capfd.readouterr()
    assert 'Axom information' in captured.out


def test_getversion_and_gitsha():
    v = pycore.getVersion()
    g = pycore.gitSHA()
    assert isinstance(v, str)
    assert isinstance(g, str)
