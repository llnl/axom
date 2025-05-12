import my_ext

# content of test_my_ext.py
def func(x):
    return x + 1


def test_answer():
    assert func(4) == 5

def test_my_ext():
    assert my_ext.add(1,4) == 5
