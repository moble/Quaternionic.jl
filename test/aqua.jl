@testitem "Aqua quality assurance" tags=[:unit] begin
    using Aqua
    Aqua.test_all(Quaternionic)
end
