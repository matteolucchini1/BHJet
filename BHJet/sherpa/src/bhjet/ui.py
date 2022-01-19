"""
Provide access to the BHJet model using the Sherpa UI layer.
"""

import logging

import sherpa.astro.ui

import bhjet

logger = logging.getLogger("sherpa")

# Tell the UI layer about the model
#
sherpa.astro.ui.add_model(bhjet.BHJet)
logger.info("Adding additive model: bhjet")
del logger
