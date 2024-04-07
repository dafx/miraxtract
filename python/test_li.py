import unittest
import librosa

def generate_features(aud_file):
    y, sr = librosa.load(aud_file)
    y = librosa.resample(y, orig_sr=sr, target_sr=12000)
    sr = 12000
    melspec = librosa.feature.melspectrogram(y=y, sr=sr, n_mels=64)
    return melspec, sr


class MyTestCase(unittest.TestCase):
    def test_something(self):
        # TODO: Write your test case here
        pass

    # BEGIN: Generated Test Cases
    def test_case1(self):
        # TODO: Write your test case here
        pass

    def test_case2(self):
        # TODO: Write your test case here
        pass
    # END: Generated Test Cases

if __name__ == '__main__':
    unittest.main()