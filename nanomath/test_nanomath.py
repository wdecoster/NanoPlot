import unittest
import nanomath as nm


class NanomathTest(unittest.TestCase):
    def test_ave_qual(self):
        """Test average quality calculation."""
        quals = list(range(128 + 1)) * 100
        mq = nm.ave_qual(quals, qround=True)
        self.assertEqual(mq, 14)


if __name__ == '__main__':
    unittest.main()
