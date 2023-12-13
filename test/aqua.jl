@testset verbose=true "Aqua quality assurance tests" begin
    Aqua.test_all(
        Quaternionic;
        ambiguities=false,
        stale_deps=(;ignore=[:Requires])  # Need Requires on Julia <1.9; not loaded on ≥1.9
    )
end
