#!/bin/bash

project.py --xml numi_fhc_workflow.xml --stage reweight_numi_fhc_run1_beam --clean
project.py --xml numi_fhc_workflow.xml --stage reweight_numi_fhc_run1_beam --submit

project.py --xml numi_fhc_workflow.xml --stage reweight_numi_fhc_run1_strangeness --clean
project.py --xml numi_fhc_workflow.xml --stage reweight_numi_fhc_run1_strangeness --submit
