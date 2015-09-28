# Copyright (C) 2009-2012 Christof Buchbender
r"""
Functions to help to treat upper/lower limits of measured Values
correctly when they are plottet on either of the axis of a
plot.

Limits are trated as follows:

When noted in the original measurement or in Text:

Lower Limit = '<'
Upper Limit = '>'

In a plot they are on the Y axis:

Lower Limit = 'v'
Upper Limit = '^'

On the X axis the convention is:

Lower Limit = '<'
Upper Limit = '>'

In case it is not clear if the limit is upper or lower. The string ``"^v"``
is used for both 
"""

def changeLimitXY(limits, direction='yToX'):
    r"""
    The changeLimitXY function changes values marked accrodingly to one axis
    to the other notation and also to and fro text notation.

    Parameters
    ----------

    limits: list
        A list of limits of arbitrary lenght. The symbols used are:

            * ``o``: No upper/lower limit, i.e. a clear detected value.
            * ``>``: A lower limit on the X Axis or in Text notation.
            * ``<`` :An upper limit on the X Axis or in Text notation
            * ``^``: An upper limit on the Y Axis.
              (Y Axis is default for PLots)
            * ``v``: An lower limit on the Y Axis.
    direction: string
        The direction of the conversion of the limits.
        Can be:

            * ``yToX`` or ``toText``: Conversion from limits for the
              y-Axis to corresponding limits on the x-Axis or for text
              notation.
            * ``xToY`` or ``toPlot``: Conversion from limits for the x-Axis
              or for text Notation to the correspondig limits for the y-Axis.

    Returns
    -------

    limits: list
        The converted limits.
    """
    if direction == 'yToX' or direction == 'toText':
        for i in range(len(limits)):
            if limits[i] == 'o':
                limits[i] = 'o'
            if limits[i] == '^':
                limits[i] = '>'
            if limits[i] == 'v':
                limits[i] = '<'
            if limits[i] == '^v':
                limits[i] = '<>'
    if direction == 'xToY' or direction == 'toPlot':
        for i in range(len(limits)):
            if limits[i] == 'o':
                limits[i] = 'o'
            if limits[i] == '>':
                limits[i] = '^'
            if limits[i] == '<':
                limits[i] = 'v'
            if limits[i] == '<>':
                limits[i] = '^v'
    return limits


def addSubLimits(first, second):
    r"""
    TODO: Remember the purpose of this function. Give better names and
    documentation in the future.
    """
    limits = []
    for i in range(len(first)):

        if first[i] == 'v' and second[i] == 'v':
            limits += ['v']
        if first[i] == '^' and second[i] == '^':
            limits += ['^']

        if first[i] == '^' and second[i] == 'v':
             limits += ['^v']
        if first[i] == 'v' and second[i] == '^':
             limits += ['^v']

        if first[i] == '^' and second[i] == 'o':
            limits += ['^']
        if first[i] == 'o' and second[i] == '^':
            limits += ['^']

        if first[i] == 'v' and second[i] == 'o':
            limits += ['v']
        if first[i] == 'o' and second[i] == 'v':
            limits += ['v']
        if first[i] == 'o' and second[i] == 'o':
            limits += ['o']
    return limits


def ratioLimits(numerator, denominator, mode='yaxis'):
    r"""
    Evaluation of upper/lower limits when building the ratio
    of two values with limits.

    Parameters
    ----------

    numerator: list
        List of arbitrary lenght of limits for the values in the nominator of
        the ratio.
    denominator: list
        List of arbitrary lenght of limits for the values in the denominator of
        the ratio. *MUST* be the same length as numerator, though.
    mode: string
        Maybe ``yaxis``(Default) or ``text`` to control the type of notation.
        See module description.
    """
    limits = []
    for i in range(len(numerator)):

        # In "Y-Axis" notation
        if mode == 'yaxis':
            if numerator[i] == '^v' or denominator[i] == '^v':
                limits += ['^v']

            if numerator[i] == 'v' and denominator[i] == 'v':
                limits += ['^v']
            if numerator[i] == '^' and denominator[i] == '^':
                limits += ['^v']

            if numerator[i] == '^' and denominator[i] == 'v':
                 limits += ['^']
            if numerator[i] == 'v' and denominator[i] == '^':
                 limits += ['v']

            if numerator[i] == '^' and denominator[i] == 'o':
                limits += ['^']
            if numerator[i] == 'o' and denominator[i] == '^':
                limits += ['v']

            if numerator[i] == 'v' and denominator[i] == 'o':
                limits += ['v']
            if numerator[i] == 'o' and denominator[i] == 'v':
                limits += ['^']
            if numerator[i] == 'o' and denominator[i] == 'o':
                limits += ['o']

        if mode == 'text':
            if numerator[i] == '^v' or denominator[i] == '^v':
                limits += ['^v']

            if numerator[i] == '<' and denominator[i] == '<':
                limits += ['^v']
            if numerator[i] == '>' and denominator[i] == '>':
                limits += ['^v']

            if numerator[i] == '>' and denominator[i] == '<':
                 limits += ['>']
            if numerator[i] == '<' and denominator[i] == '>':
                 limits += ['<']

            if numerator[i] == '>' and denominator[i] == 'o':
                limits += ['>']
            if numerator[i] == 'o' and denominator[i] == '>':
                limits += ['<']

            if numerator[i] == '<' and denominator[i] == 'o':
                limits += ['<']
            if numerator[i] == 'o' and denominator[i] == '<':
                limits += ['>']

            if numerator[i] == 'o' and denominator[i] == 'o':
                limits += ['o']
    return limits


def limCheck(xAxisLimits, yAxisLimits):
    r""" This function determines which upper/lower limits have to be set on
    the x and y Axis. This is helpful when working with matplotlib.

    TODO: Extend documentation.
    It uses False, True arguments in a list of
    four entries the entries correspond to: limitChecker=
    [xAxisLow, xAxisUp, yAxisLow, yAxisUp] by default the values are
    set to False.
    """
    limitChecker=[False, False, False, False]
    if yAxisLimits == 'v':
        limitChecker[2] = True
    if yAxisLimits == '^':
        limitChecker[3] = True
    if yAxisLimits == '^v':
        limitChecker[2] = True
        limitChecker[3] = True

    if xAxisLimits == '<':
        limitChecker[0] = True
    if xAxisLimits == '>':
        limitChecker[1] = True
    if xAxisLimits == '<>':
        limitChecker[0] = True
        limitChecker[1] = True
    return limitChecker
