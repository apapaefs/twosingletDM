import tempfile
import unittest
from pathlib import Path

from scan_output import write_valid_point


class TestScanOutput(unittest.TestCase):
    def test_write_valid_point_adds_header_and_mg5_xsecs(self):
        point_info = {
            "M2": 300.0,
            "M3": 500.0,
            "vs": 750.0,
            "vx": 0.0,
            "a12": 0.1,
            "a13": 0.0,
            "a23": 0.0,
            "lX": 0.2,
            "lPhiX": 0.3,
            "lSX": 0.4,
            "w1": 0.004,
            "w2": 1.1,
            "w3": 2.2,
            "k1": 0.99,
            "k2": 0.1,
            "k3": 0.0,
            "K111": 1.0,
            "K112": 2.0,
            "K113": 3.0,
            "K123": 4.0,
            "K122": 5.0,
            "K1111": 6.0,
            "K1112": 7.0,
            "K1113": 8.0,
            "K133": 9.0,
            "evo": True,
            "thc": True,
            "hb": True,
            "hs": False,
            "ewpo": True,
            "wmass": True,
            "dm": True,
            "dm_mdm": 1000.0,
            "dm_omega": 0.049,
            "dm_relic_upper_limit": 0.1224,
            "dm_dir_det": 1.218e-9,
            "dm_dir_det_limit": 2.05907e-9,
            "dm_lux_base_limit": 8.236e-10,
            "dm_relic_excluded": False,
            "dm_direct_detection_excluded": False,
        }

        with tempfile.TemporaryDirectory() as tmpdir:
            outfile = Path(tmpdir) / "points.dat"
            write_valid_point(outfile, point_info, {"hh": 1.25, "hhh": 0.031})

            lines = outfile.read_text(encoding="ascii").splitlines()

        header = lines[0].split("\t")
        row = lines[1].split("\t")
        self.assertEqual(len(header), len(row))
        self.assertEqual(header[:4], ["M2", "M3", "vs", "vx"])
        self.assertIn("mg5_xsec_hh_pb", header)
        self.assertIn("mg5_xsec_hhh_pb", header)
        self.assertIn("dm_omega", header)
        self.assertIn("dm_dir_det_limit", header)
        self.assertEqual(row[header.index("M2")], "300.0")
        self.assertEqual(row[header.index("hs")], "False")
        self.assertEqual(row[header.index("dm_omega")], "0.049")
        self.assertEqual(row[header.index("dm_dir_det_limit")], "2.05907e-09")
        self.assertEqual(row[header.index("mg5_xsec_hh_pb")], "1.25")
        self.assertEqual(row[header.index("mg5_xsec_hhh_pb")], "0.031")

    def test_write_valid_point_does_not_repeat_header(self):
        point_info = {"M2": 300.0, "M3": 500.0, "vs": 750.0, "vx": 0.0}

        with tempfile.TemporaryDirectory() as tmpdir:
            outfile = Path(tmpdir) / "points.dat"
            write_valid_point(outfile, point_info, {})
            write_valid_point(outfile, point_info, {})

            lines = outfile.read_text(encoding="ascii").splitlines()

        self.assertEqual(len(lines), 3)


if __name__ == "__main__":
    unittest.main()
