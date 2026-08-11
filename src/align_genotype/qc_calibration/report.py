"""Human-readable formatting for the calibration commands."""


def fmt_measure(value: float, unit: str) -> str:
    """Render a metric value at a precision that reads sensibly for its unit."""
    if unit == 'frac':
        return f'{value:.3f}'
    if unit in ('x', '%'):
        return f'{value:.1f}'
    return f'{value:.2f}'
