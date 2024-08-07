from __future__ import division, absolute_import, print_function

import pytest

from os.path import dirname, join, abspath

from cobrame.io.json import (load_reduced_json_me_model, save_json_me_model,
                             load_json_me_model)
from cobrame.util.massbalance import check_me_model_mass_balance
from ecolime.build_me_model import return_me_model
from ecolime.util.me_model_comparison import find_me_model_difference

current_dir = dirname(abspath(__file__))
models_dir = current_dir.split('tests')[0] + 'me_models/'

del dirname, abspath

test_model = return_me_model()
save_json_me_model(test_model, 'test_json_dump.json')
json_model = load_json_me_model('test_json_dump.json')


def test_model_benchmark():
    benchmark_model = \
        load_reduced_json_me_model(join(models_dir, 'iJL1678b_reduced.json'))
    difference = find_me_model_difference(benchmark_model, test_model, 1e-10)
    print('-----------------------Difference----------------------')
    print(difference)
    assert (len(difference) == 0)


def test_full_json_dumping():
    difference = find_me_model_difference(test_model, json_model, 1e-10)
    print(difference)
    assert (len(difference) == 0)


def test_mass_balance():
    not_mass_balanced = check_me_model_mass_balance(test_model)
    not_mass_balanced.pop('23bpg_generation_FWD_CPLX_dummy')
    not_mass_balanced.pop('GLUTRR_FWD_CPLX0-3741')
    not_mass_balanced.pop('dpm_import_FWD_CPLX_dummy')
    not_mass_balanced.pop('tl_generation_FWD_CPLX_dummy')
    print(not_mass_balanced)
    assert (len(not_mass_balanced) == 0)


def test_json_mass_balance():
    not_mass_balanced = check_me_model_mass_balance(json_model)
    not_mass_balanced.pop('23bpg_generation_FWD_CPLX_dummy')
    not_mass_balanced.pop('GLUTRR_FWD_CPLX0-3741')
    not_mass_balanced.pop('dpm_import_FWD_CPLX_dummy')
    not_mass_balanced.pop('tl_generation_FWD_CPLX_dummy')
    print(not_mass_balanced)
    assert (len(not_mass_balanced) == 0)


if __name__ == '__main__':
    test_model_benchmark()
