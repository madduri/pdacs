#!/usr/bin/env python


def ds_inputOptions(inputmodel,inputSize,inputRealization):
    """List available input Models as tuples of (displayName,value)"""
    print inputmodel + " " + inputSize + " " + inputRealization
    models = [["g","g","g"],["h","h","h"],["i","i","i"]]
    return models
#      <param name="inputmodel" type="select" size="100" label="Input Model" area="True" dynamic_options="ds_inputmodelOptions()"/>
#      <param name="inputSize" type="select" label="Box Size" selected="false" dynamic_options="ds_inputsizeOptions()"/>
#      <param name="inputRealization" type="select" label="Realization based on different seed values" selected="false" dynamic_options="ds_inputrealizationOptions()"/>
#      <param name="input" type="select" size="300" label="Input Snapshot" area="True" help="Choosing Input Snapshot files" selected="false" dynamic_options="ds_inputOptions()"/>

