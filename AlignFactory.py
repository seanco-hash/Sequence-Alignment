from GlobalAlign import GlobalAlign
from OverlapAlign import OverlapAlign
from LocalAlign import LocalAlign


def getAlign(s, t, align_type, scores):
    """
    This factory creates alignments from the following types: global,
    local and overlap.
    :param s: sequence 1
    :param t: sequence 2
    :param align_type: wanted type
    :param scores: scores matrix which the alignment works by.
    :return: alignment instance
    """
    if align_type == "overlap":
        return OverlapAlign(s, t, align_type, scores)
    if align_type == "global":
        return GlobalAlign(s, t, align_type, scores)
    return LocalAlign(s, t, align_type, scores)
