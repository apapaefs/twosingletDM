if __package__:
    from .flavour_observables import *
else:
    from flavour_observables import *


if __name__ == "__main__":

    # Example parameter point
    m_hi_test = 5.0
    R_ih_test = 0.01
    BR_mumu_test = 0.10
    BR_tautau_test = 0.80

    production = BRUpsgammahi(
        m_hi=m_hi_test,
        R_ih=R_ih_test,
    )

    prediction_mumu = BRUpsgammahiFinalState(
        m_hi=m_hi_test,
        R_ih=R_ih_test,
        BR_hi_to_final=BR_mumu_test,
    )

    prediction_tautau = BRUpsgammahiFinalState(
        m_hi=m_hi_test,
        R_ih=R_ih_test,
        BR_hi_to_final=BR_tautau_test,
    )

    limit_mumu = min(
        belle_mumu_limit(m_hi_test),
        babar_mumu_limit(m_hi_test),
    )

    limit_tautau = min(
        belle_tautau_limit(m_hi_test),
        babar_tautau_limit(m_hi_test),
    )

    result = CheckUpsLeptonBounds(
        m_hi=m_hi_test,
        R_ih=R_ih_test,
        BR_hi_to_mumu=BR_mumu_test,
        BR_hi_to_tautau=BR_tautau_test,
    )

    print(f"Scalar mass                   = {m_hi_test:.3f} GeV")
    print(f"R_ih                          = {R_ih_test}")
    print(f"BR(h_i -> mu+ mu-)            = {BR_mumu_test:.6e}")
    print(f"BR(h_i -> tau+ tau-)          = {BR_tautau_test:.6e}")
    print()
    print(f"BR(Upsilon -> gamma h_i)      = {production:.6e}")
    print()
    print("--- Muon channel ---")
    print(f"Predicted product BR          = {prediction_mumu:.6e}")
    print(f"Belle limit                   = {belle_mumu_limit(m_hi_test):.6e}")
    print(f"BaBar limit                   = {babar_mumu_limit(m_hi_test):.6e}")
    print(f"Strongest limit               = {limit_mumu:.6e}")
    print()
    print("--- Tau channel ---")
    print(f"Predicted product BR          = {prediction_tautau:.6e}")
    print(f"Belle limit                   = {belle_tautau_limit(m_hi_test):.6e}")
    print(f"BaBar limit                   = {babar_tautau_limit(m_hi_test):.6e}")
    print(f"Strongest limit               = {limit_tautau:.6e}")
    print()
    print(f"Constraint result             = {result}")
    print("Point allowed" if result == 1 else "Point excluded")
