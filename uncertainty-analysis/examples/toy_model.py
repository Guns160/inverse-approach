"""Synthetic algebraic model for exercising the two separated ensembles."""


def evaluate(parameters):
    if parameters["scale"] <= 0 or parameters["Y2"] <= parameters["Y1"]:
        raise ValueError("invalid synthetic model domain")
    return {
        "sum_response": (parameters["Q2"] + parameters["condensate"]
                         + 0.01 * (parameters["Y2"] - parameters["Y1"])
                         + 0.1 * parameters["X2"]),
        "window_width": parameters["Y2"] - parameters["Y1"],
        "scaled_product": (parameters["Q2"] * parameters["condensate"]
                            * parameters["m0"] / parameters["scale"]),
    }
