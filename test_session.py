import pytest
from session import add_numbers

# test positive
def test_add_positive():
    assert add_numbers(1,2) == 3

# test add 0
def test_add_zero():
    assert add_numbers(1,0) == 1

# test negative
def test_add_negative():
    assert add_numbers(4,-100) == -96

# type error 
def test_add_string_expect_exception():
    with pytest.raises(TypeError):
        add_numbers(4, 'Anything')
