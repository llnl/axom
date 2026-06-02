import pycore

def test_about_returns_string():
    s = pycore.about()
    assert isinstance(s, str)
    assert 'Axom information' in s


def test_getversion_and_gitsha():
    v = pycore.getVersion()
    g = pycore.gitSHA()
    assert isinstance(v, str)
    assert isinstance(g, str)
