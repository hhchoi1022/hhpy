


class Helper:
    
    def __init__(self):
        from Research.helper import PhotometryHelper
        from Research.helper import SpectroscopyHelper
        from Research.helper import AnalysisHelper
        self.photometry = PhotometryHelper()
        self.spectroscopy = SpectroscopyHelper()
        self.analyze = AnalysisHelper()
        pass
    