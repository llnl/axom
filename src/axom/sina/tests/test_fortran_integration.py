import argparse
import io
import json
import os
import subprocess
import unittest
import config

if (config.AXOM_USE_HDF5):
    import h5py


def parse_args():
    """Helper function to obtain the binary directory path of Axom from CLI"""
    parser = argparse.ArgumentParser(description="Unit test arguments")
    parser.add_argument("-bd", "--binary-dir", type=str,
                        help="Path to the binary directory for Axom")
    # Add other arguments as needed
    return parser.parse_args()

# JSON Tests: Will always run
class TestFortranExampleIntegrationJSON(unittest.TestCase):
    
    @classmethod
    def setUpClass(cls):
        """
        Obtain the binary directory from the CLI and compile the sina fortran
        example needed for these tests if necessary.
        """
        cwd = os.getcwd()

        args = parse_args()
        cls.binary_dir = args.binary_dir
        if cls.binary_dir is None:
            # Move up three levels and resolve an absolute path
            cls.binary_dir = os.path.abspath(os.path.join(cwd, "../"))

        os.chdir(cls.binary_dir)

        if not os.path.exists(os.path.join(cls.binary_dir, "examples/sina_fortran_ex")):
            subprocess.run(["make", "sina_fortran_ex"])

        os.chdir(cwd)

    def setUp(self):
        """ Invoke example Fortran application to dump a sina file """
        subprocess.run([os.path.join(self.binary_dir, "examples/sina_fortran_ex")])
        self.dump_file = "sina_dump.json"
        
    def tearDown(self):
        """ Clean up output directory after each test. """
        os.remove(self.dump_file)
    
    def test_file_validity(self):
        """ Make sure the files we're importing follow the Sina schema. """
        try:
            import jsonschema
            schema_file = os.path.join(self.binary_dir, "tests", "sina_schema.json")
            
            with io.open(schema_file, "r", encoding="utf-8") as schema:
                schema = json.load(schema)
                with io.open(self.dump_file, "r", encoding="utf-8") as loaded_test:
                    file_json = json.load(loaded_test)
                    jsonschema.validate(file_json, schema)
        except ModuleNotFoundError:
            print("jsonschema module not found. Skipping test_file_validity.")
            pass
                    
    def test_validate_contents_of_record(self):
        """ Ensure that the  record written out matches what we expect """
        with open(self.dump_file, "r", encoding="utf-8") as loaded_test:
            rec = json.load(loaded_test)

        record = rec["records"][0]
        
        # Test the metadata in the record
        self.assertEqual("my_rec_id", record["id"])
        self.assertEqual("my_type", record["type"])
        
        # Test the files
        self.assertEqual(list(record["files"].keys()), ["/path/to/my/file/my_other_file.txt", "/path/to/my/file/my_file.txt"])
        self.assertEqual(record["files"]["/path/to/my/file/my_other_file.txt"]["mimetype"], "png")
        self.assertEqual(record["files"]["/path/to/my/file/my_file.txt"]["mimetype"], "txt")
        
        # Test the signed variants 
        self.assertEqual("A", record["data"]["char"]["value"])
        self.assertEqual(10, record["data"]["int"]["value"])
        self.assertEqual(0, record["data"]["logical"]["value"])
        self.assertEqual(1000000000.0, record["data"]["long"]["value"])
        self.assertEqual(1.23456704616547, record["data"]["real"]["value"])
        self.assertEqual(0.810000002384186, record["data"]["double"]["value"])
        
        # Test the unsigned variants
        self.assertEqual("A", record["data"]["u_char"]["value"])
        self.assertEqual("kg", record["data"]["u_char"]["units"])
        self.assertEqual(10, record["data"]["u_int"]["value"])
        self.assertEqual("kg", record["data"]["u_int"]["units"])
        self.assertEqual(1.0, record["data"]["u_logical"]["value"])
        self.assertEqual("kg", record["data"]["u_logical"]["units"])
        self.assertEqual(1000000000.0, record["data"]["u_long"]["value"])
        self.assertEqual("kg", record["data"]["u_long"]["units"])
        self.assertEqual(1.23456704616547, record["data"]["u_real"]["value"])
        self.assertEqual("kg", record["data"]["u_real"]["units"])
        self.assertEqual(0.810000002384186, record["data"]["u_double"]["value"])
        self.assertEqual("kg", record["data"]["u_double"]["units"])
        self.assertEqual(0.810000002384186, record["data"]["u_double_w_tag"]["value"])
        self.assertEqual("kg", record["data"]["u_double_w_tag"]["units"])
        self.assertEqual(["new_fancy_tag"], record["data"]["u_double_w_tag"]["tags"])
        
        # Test the curves
        nums = range(1, 21)
        real_arr = [i for i in nums]
        double_arr = [i*2 for i in nums]
        int_arr = [i*3 for i in nums]
        long_arr = [i*4 for i in nums]
        curveset = "my_curveset"
        for kind, loc in (("indep", "independent"), ("dep", "dependent")):
            for val_type, target in (("real", real_arr), ("double", double_arr), ("int", int_arr), ("long", long_arr)):
                name = "my_{}_curve_{}".format(kind, val_type)
                self.assertEqual(target, record["curve_sets"][curveset][loc][name]["value"])
        double_2_name = "my_dep_curve_double_2"
        self.assertEqual(double_arr, record["curve_sets"][curveset]["dependent"][double_2_name]["value"])


#HDF5 Test
@unittest.skipUnless(config.AXOM_USE_HDF5, "Requires h5py for HDF5-dependent tests")
class TestFortranExampleIntegrationHDF5(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        """
        Obtain the binary directory from the CLI and compile the sina fortran
        example needed for these tests if necessary.
        """
        cwd = os.getcwd()

        args = parse_args()
        cls.binary_dir = args.binary_dir
        if cls.binary_dir is None:
            # Move up three levels and resolve an absolute path
            cls.binary_dir = os.path.abspath(os.path.join(cwd, "../"))

        os.chdir(cls.binary_dir)

        if not os.path.exists(os.path.join(cls.binary_dir, "examples/sina_fortran_ex")):
            subprocess.run(["make", "sina_fortran_ex"])

        os.chdir(cwd)

    def setUp(self):
        """ Invoke example Fortran application to dump a sina file """
        subprocess.run([os.path.join(self.binary_dir, "examples/sina_fortran_ex")])
        self.dump_file = "sina_dump.hdf5"

    def tearDown(self):
        os.remove(self.dump_file)

    def extract_hdf5_value(self, value):
        # If the value is an h5py.Dataset, retrieve its underlying data.
        if isinstance(value, h5py.Dataset):
            value = value[()]
            
        # If the value is a bytes instance, decode it.
        if isinstance(value, bytes):
            return value.decode("utf-8").strip("\0").strip()
        
        # If the value is a list or tuple of bytes, join them into a single bytes object and decode.
        if isinstance(value, (list, tuple)) and value and all(isinstance(item, bytes) for item in value):
            joined = b"".join(value)
            return joined.decode("utf-8").strip("\0").strip()
        
        # If the value has a 'tolist' method (e.g., from an array-like object), convert it and process.
        if hasattr(value, "tolist"):
            converted = value.tolist()
            if isinstance(converted, (list, tuple)) and converted and all(isinstance(item, bytes) for item in converted):
                joined = b"".join(converted)
                return joined.decode("utf-8").strip("\0").strip()
            if isinstance(converted, list) and len(converted) == 1:
                return converted[0]
            return converted

        # Otherwise, return the value as-is.
        return value

    def test_validate_contents_of_record(self):
        with h5py.File(self.dump_file, "r") as f:
            record = f["records"]["0"]

            # Validate metadata
            self.assertEqual("my_rec_id", self.extract_hdf5_value(record["id"]))
            self.assertEqual("my_type", self.extract_hdf5_value(record["type"]))

            # Validate Files
            files_group = record["files"]
            expected_file_keys = [
                "__SINA_SLASHREPLACE__path__SINA_SLASHREPLACE__to__SINA_SLASHREPLACE__my__SINA_SLASHREPLACE__file__SINA_SLASHREPLACE__my_other_file.txt",
                "__SINA_SLASHREPLACE__path__SINA_SLASHREPLACE__to__SINA_SLASHREPLACE__my__SINA_SLASHREPLACE__file__SINA_SLASHREPLACE__my_file.txt"
            ]
            self.assertEqual(sorted(list(files_group.keys())), sorted(expected_file_keys))
            self.assertEqual(self.extract_hdf5_value(
                files_group["__SINA_SLASHREPLACE__path__SINA_SLASHREPLACE__to__SINA_SLASHREPLACE__my__SINA_SLASHREPLACE__file__SINA_SLASHREPLACE__my_other_file.txt"]["mimetype"]),
                "png")
            self.assertEqual(self.extract_hdf5_value(
                files_group["__SINA_SLASHREPLACE__path__SINA_SLASHREPLACE__to__SINA_SLASHREPLACE__my__SINA_SLASHREPLACE__file__SINA_SLASHREPLACE__my_file.txt"]["mimetype"]),
                "txt")

            # Validate Data
            data_group = record["data"]
            self.assertEqual(self.extract_hdf5_value(data_group["char"]["value"]), "A")
            self.assertEqual(self.extract_hdf5_value(data_group["int"]["value"]), 10)
            self.assertEqual(self.extract_hdf5_value(data_group["logical"]["value"]), 0)
            self.assertEqual(self.extract_hdf5_value(data_group["long"]["value"]), 1000000000.0)
            self.assertAlmostEqual(self.extract_hdf5_value(data_group["real"]["value"]), 1.23456704616547)
            self.assertAlmostEqual(self.extract_hdf5_value(data_group["double"]["value"]), 0.810000002384186)

            self.assertEqual(self.extract_hdf5_value(data_group["u_char"]["value"]), "A")
            self.assertEqual(self.extract_hdf5_value(data_group["u_char"]["units"]), "kg")
            self.assertEqual(self.extract_hdf5_value(data_group["u_int"]["value"]), 10)
            self.assertEqual(self.extract_hdf5_value(data_group["u_int"]["units"]), "kg")
            self.assertEqual(self.extract_hdf5_value(data_group["u_logical"]["value"]), 1.0)
            self.assertEqual(self.extract_hdf5_value(data_group["u_logical"]["units"]), "kg")
            self.assertEqual(self.extract_hdf5_value(data_group["u_long"]["value"]), 1000000000.0)
            self.assertEqual(self.extract_hdf5_value(data_group["u_long"]["units"]), "kg")
            self.assertAlmostEqual(self.extract_hdf5_value(data_group["u_real"]["value"]), 1.23456704616547)
            self.assertEqual(self.extract_hdf5_value(data_group["u_real"]["units"]), "kg")
            self.assertAlmostEqual(self.extract_hdf5_value(data_group["u_double"]["value"]), 0.810000002384186)
            self.assertEqual(self.extract_hdf5_value(data_group["u_double"]["units"]), "kg")
            self.assertAlmostEqual(self.extract_hdf5_value(data_group["u_double_w_tag"]["value"]), 0.810000002384186)
            self.assertEqual(self.extract_hdf5_value(data_group["u_double_w_tag"]["units"]), "kg")

            tags_value = data_group["u_double_w_tag"]["tags"]
            if isinstance(tags_value, h5py.Group):
                if len(tags_value) == 1:
                    inner_dataset = list(tags_value.values())[0]
                    tags_value = self.extract_hdf5_value(inner_dataset)
            if isinstance(tags_value, str):
                tags_value = [tag.strip() for tag in tags_value.split(',')]
            self.assertEqual(tags_value, ["new_fancy_tag"])

            # Validate Curves
            curveset_group = record["curve_sets"]["my_curveset"]
            independent_group = curveset_group["independent"]
            dependent_group = curveset_group["dependent"]

            nums = list(range(1, 21))
            real_arr   = [i for i in nums]
            double_arr = [i * 2 for i in nums]
            int_arr    = [i * 3 for i in nums]
            long_arr   = [i * 4 for i in nums]

            for kind, grp in (("indep", independent_group), ("dep", dependent_group)):
                for val_type, target in (("real", real_arr),
                                         ("double", double_arr),
                                         ("int", int_arr),
                                         ("long", long_arr)):
                    curve_name = f"my_{kind}_curve_{val_type}"
                    self.assertIn(curve_name, grp)
                    curve_val = self.extract_hdf5_value(grp[curve_name]["value"])
                    self.assertEqual(curve_val, target)

            double_2_name = "my_dep_curve_double_2"
            self.assertIn(double_2_name, dependent_group)
            curve_double_2 = self.extract_hdf5_value(dependent_group[double_2_name]["value"])
            self.assertEqual(curve_double_2, double_arr)


if __name__ == "__main__":
    # Doing the below instead of unittest.main() so that we can print to stdout
    suite = unittest.TestSuite()
    suite.addTests(unittest.defaultTestLoader.loadTestsFromTestCase(TestFortranExampleIntegrationJSON))
    if config.AXOM_USE_HDF5:
        suite.addTests(unittest.defaultTestLoader.loadTestsFromTestCase(TestFortranExampleIntegrationHDF5))
    runner = unittest.TextTestRunner(buffer=False)
    runner.run(suite)
