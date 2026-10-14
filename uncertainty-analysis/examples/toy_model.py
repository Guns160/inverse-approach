"""Synthetic algebraic example only: no physical OPE or numerical extraction."""

def evaluate(parameters):
    q2 = parameters["Q2"]
    condensate = parameters["C_example"]
    scale = parameters["scale"]
    if scale <= 0:
        raise ValueError("scale must be positive")
    return {
        "sum_response": q2 + condensate,
        "difference_response": q2 - condensate,
        "scaled_product": q2 * condensate / scale,
    }
