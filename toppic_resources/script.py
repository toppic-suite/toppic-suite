def fit_vogit(apex, yhat):
    import numpy
    from lmfit.models import SkewedVoigtModel
    xxx = numpy.linspace(0, len(yhat) - 1, len(yhat))
    v_list_filtered = numpy.array([val for val in yhat])
    model = SkewedVoigtModel()
    params = model.guess(v_list_filtered, x=xxx)
    params['center'].set(value=apex)
    result = model.fit(v_list_filtered, params, x=xxx)
    fitted_data = result.best_fit.tolist()
    return fitted_data
