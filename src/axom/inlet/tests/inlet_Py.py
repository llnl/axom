from pathlib import Path

import pyinlet


def make_inlet(reader, source: str):
    reader.parseString(source)
    return pyinlet.Inlet(reader)


def test_get_nested_values_added_using_container():
    inlet = make_inlet(
        pyinlet.YAMLReader(),
        "foo:\n  bar: yet another string\n  so: 3.5\n  re: 9\n  mi: 1\n",
    )

    container = inlet.addStruct("foo", "A container called foo")
    container.required(True)

    container.addString("bar", "bar's description").required(True)
    container.addDouble("so", "so's description")
    container.addInt("re", "re's description")
    container.addInt("mi", "mi's description")

    assert str(inlet["foo/bar"]) == "yet another string"
    assert int(inlet["foo/mi"]) == 1
    assert float(inlet["foo/so"]) == 3.5
    assert int(inlet["foo/re"]) == 9


def test_get_nested_ints_through_container():
    inlet = make_inlet(pyinlet.YAMLReader(), "foo:\n  bar: 1\n  baz: 2\n")
    foo = inlet.addStruct("foo", "A container called foo")
    foo.addInt("bar", "bar's description")
    foo.addInt("baz", "baz's description")

    assert int(foo["bar"]) == 1
    assert int(foo["baz"]) == 2


def test_simple_struct_verify_pass():
    inlet = make_inlet(pyinlet.YAMLReader(), "foo:\n  bar: 1\n  baz: 2\n")
    foo = inlet.addStruct("foo", "A struct called foo")
    foo.addInt("bar", "bar's description").required(True)
    foo.addInt("baz", "baz's description").required(True)

    assert inlet.verify()


def test_simple_struct_verify_fail():
    inlet = make_inlet(pyinlet.YAMLReader(), "foo:\n  baz: 2\n")
    foo = inlet.addStruct("foo", "A struct called foo")
    foo.addInt("bar", "bar's description").required(True)
    foo.addInt("baz", "baz's description").required(True)

    assert not inlet.verify()


def test_default_scalar_value():
    reader = pyinlet.YAMLReader()
    reader.parseString(" ")
    inlet = pyinlet.Inlet(reader)

    foo = inlet.addStruct("foo", "A container called foo")
    foo.addInt("bar", "bar's description").defaultValue(2)
    foo.addDouble("baz", "baz's description").defaultValue(3.14)
    foo.addString("quux", "quux's description").defaultValue("default string")

    assert int(inlet["foo/bar"]) == 2
    assert float(inlet["foo/baz"]) == 3.14
    assert str(inlet["foo/quux"]) == "default string"
    assert inlet.contains("foo")
    assert inlet.verify()


def test_parse_file(tmp_path: Path):
    input_file = tmp_path / "input.yaml"
    input_file.write_text("foo:\n  bar: 1\n", encoding="utf-8")

    reader = pyinlet.YAMLReader()
    assert reader.parseFile(str(input_file))
    inlet = pyinlet.Inlet(reader)
    foo = inlet.addStruct("foo", "A container called foo")
    foo.addInt("bar", "bar's description")

    assert inlet.contains("foo")
    assert int(foo["bar"]) == 1


def test_lua_reader_available_if_enabled():
    if not hasattr(pyinlet, "LuaReader"):
        return

    inlet = make_inlet(pyinlet.LuaReader(), "foo = { bar = 1; baz = 2 }")
    foo = inlet.addStruct("foo", "A container called foo")
    foo.addInt("bar", "bar's description")
    foo.addInt("baz", "baz's description")

    assert int(foo["bar"]) == 1
    assert int(foo["baz"]) == 2
