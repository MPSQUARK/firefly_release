import sys
sys.path.append('.')
from src.setup import Setup
from src.config import Config

config = Config().create("config/default.yaml")
setup = Setup(config)

def test_mask_emissionlines():
    setup.mask_emissionlines(config.emlines)
    pass
