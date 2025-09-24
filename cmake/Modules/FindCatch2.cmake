include(FetchContent)

# Catch2 v2.x branch (i.e., version compatible with ApprovalTests)
FetchContent_Declare(
        Catch2
        GIT_REPOSITORY https://github.com/catchorg/Catch2.git
        GIT_TAG ee1450f268dfd5c13aa8670ba97e93cabaf2e15d
)
FetchContent_MakeAvailable(Catch2)