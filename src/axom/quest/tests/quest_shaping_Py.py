from pathlib import Path

import pyquest


def find_repo_root() -> Path:
    for candidate in [Path.cwd(), *Path.cwd().parents]:
        if (candidate / "src/examples/shaping_tutorial/lesson_04/circle_input.lua").exists():
            return candidate
    raise RuntimeError("Could not locate the Axom source tree for the quest shaping test.")


def test_run_shaping(tmp_path: Path):
    repo_root = find_repo_root()
    mesh_file = repo_root / "src/examples/shaping_tutorial/lesson_04/circle_input.lua"
    klee_file = repo_root / "src/examples/shaping_tutorial/lesson_04/circles.yaml"

    result = pyquest.runShaping(str(mesh_file), str(klee_file), str(tmp_path / "circle"), False)
    result.save()

    assert result.getBlueprintGroup() is not None
    assert result.getDataStore() is not None
