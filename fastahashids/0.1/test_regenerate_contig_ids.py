import unittest

script = __import__("regenerate_contig_ids")


class TestSum(unittest.TestCase):
  def test_megahit_header_k141(self):
    new_id = script.extract_unique_element_from_contigid("k141_5432123")
    self.assertEqual(new_id, "5432123")

  def test_megahit_header_not_k141(self):
    new_id = script.extract_unique_element_from_contigid("k31_34534")
    self.assertEqual(new_id, "k31_34534")

  def test_megahit_header_not_k141_2(self):
    new_id = script.extract_unique_element_from_contigid("k121_1222111")
    self.assertEqual(new_id, "k121_1222111")

  def test_spades(self):
    new_id = script.extract_unique_element_from_contigid("NODE_1232_length_891644_cov_69.428772")
    self.assertEqual(new_id, "1232")

  def test_spades_2(self):
    new_id = script.extract_unique_element_from_contigid("NODE_20532_length_1001_cov_1.171247")
    self.assertEqual(new_id, "20532")

  def test_some_other_assembler_id(self):
    try:
      script.extract_unique_element_from_contigid("OTHER_ASSEMBLER_123445")
      self.assertTrue(False)
    except NotImplementedError:
      self.assertTrue(True)


if __name__ == '__main__':
  unittest.main()
