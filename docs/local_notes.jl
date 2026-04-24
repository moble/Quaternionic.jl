# Conditionally include private local notes when building docs.  The directory
# `docs/src/local_notes` should be a symlink into a private notes git repo.  When absent
# (e.g., in CI), `notes_pages` and `notes_remotes` are empty and the build proceeds
# normally.

function local_notes()
    notes_src = joinpath(@__DIR__, "src", "local_notes")

    if isdir(notes_src)
        files = filter(f -> endswith(f, ".md"), readdir(notes_src; sort=true))
        notes_pages = isempty(files) ? [] : ["Local Notes" => map(f -> "local_notes/$f", files)]
        try
            notes_root = readchomp(`git -C $(realpath(notes_src)) rev-parse --show-toplevel`)
            notes_remote_url = readchomp(`git -C $notes_root remote get-url origin`)
            # Parse both SSH (git@github.com:user/repo.git) and HTTPS (https://github.com/user/repo.git)
            notes_slug = replace(notes_remote_url, r"^.*github\.com[:/]" => "", r"\.git$" => "")
            notes_remotes = Dict(notes_root => Documenter.Remotes.GitHub(notes_slug))
            return (notes_pages, notes_remotes)
        catch err
            @warn "Skipping local notes git remote configuration; docs build will continue without remotes for local notes." exception=(err, catch_backtrace()) notes_src=notes_src
            return (notes_pages, Dict())
        end
    else
        return ([], Dict())
    end
end
