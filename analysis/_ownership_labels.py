"""Presentation names for registration records; analytical codes stay unchanged."""

REGISTRATION_DISPLAY_NAMES = {"HKG": "Hong Kong, China"}


def normalize_registration_names(table):
    """Return a copy with canonical owner-registration display names."""
    result = table.copy()
    if {"investor_iso3", "investor_country"}.issubset(result.columns):
        for code, name in REGISTRATION_DISPLAY_NAMES.items():
            result.loc[result["investor_iso3"].eq(code), "investor_country"] = name
    return result
