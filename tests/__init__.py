# -*- coding: utf-8 -*-

import os

enptRepo_rootpath = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))

config_for_testing = dict(
    path_l1b_enmap_image=os.path.join(enptRepo_rootpath, 'tests', 'data', 'EnMAP_Level_1B', 'AlpineTest1_CWV2_SM0.zip')
)
