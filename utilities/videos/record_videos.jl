function record_videos(video_name, ego_polys, obs_polys, T, dt, R_cost, Q_cost, u1_max, u2_max, u3_max, xlims, ylims, framerate)

    nonsmooth_prob = PolyPlanning.setup_nonsmooth(
        ego_polys,
        obs_polys;
        T,
        dt,
        R_cost,
        Q_cost,
        u1_max,
        u2_max,
        u3_max,
        n_sd_slots=2
    )

    our_sol = PolyPlanning.solve_nonsmooth(nonsmooth_prob, x0; is_displaying=false, sleep_duration=0.3, is_recording=true)

    (fig, update_fig, ax) = PolyPlanning.visualize_nonsmooth(x0, T, ego_polys, obs_polys; is_displaying=true)
    PolyPlanning.xlims!(ax, xlims[1], xlims[2])
    PolyPlanning.ylims!(ax, ylims[1], ylims[2])
    PolyPlanning.hidedecorations!(ax)
    hist_len = length(our_sol.θ_history)
    frames = 2:hist_len+1*framerate

    PolyPlanning.record(fig, "videos/ours_$video_name.mp4", frames;
        framerate) do frame
        #@infiltrate
        if frame <= hist_len
            update_fig(our_sol.θ_history[frame])
        end
    end

    sep_prob = PolyPlanning.setup_sep_planes(
        ego_polys,
        obs_polys;
        T,
        dt,
        Rf=R_cost,
        Qf=Q_cost,
        u1_max,
        u2_max,
        u3_max
    )

    sep_sol = PolyPlanning.solve_prob_sep_planes(sep_prob, x0; is_displaying=false, sleep_duration=0.1, is_recording=true)

    (fig, update_fig, ax) = PolyPlanning.visualize_sep_planes(x0, T, ego_polys, obs_polys; n_per_col=2, is_displaying=false)
    PolyPlanning.xlims!(ax, xlims[1], xlims[2])
    PolyPlanning.ylims!(ax, ylims[1], ylims[2])
    PolyPlanning.hidedecorations!(ax)
    hist_len = length(sep_sol.θ_history)
    frames = 2:hist_len+1*framerate

    PolyPlanning.record(fig, "videos/sep_$video_name.mp4", frames;
        framerate) do frame
        if frame <= hist_len
            update_fig(sep_sol.θ_history[frame])
        end
    end

    dcol_prob = PolyPlanning.setup_dcol(
        ego_polys;
        T,
        dt,
        Rf=R_cost,
        Qf=Q_cost,
        u1_max,
        u2_max,
        u3_max,
        n_obs=length(obs_polys)
    )

    dcol_sol = PolyPlanning.solve_dcol(dcol_prob, x0, obs_polys; is_displaying=false, sleep_duration=0.1, is_recording=true)

    (fig, update_fig, ax) = PolyPlanning.visualize_dcol(x0, T, ego_polys, obs_polys; is_displaying=false)
    PolyPlanning.xlims!(ax, xlims[1], xlims[2])
    PolyPlanning.ylims!(ax, ylims[1], ylims[2])
    PolyPlanning.hidedecorations!(ax)
    hist_len = length(dcol_sol.θ_history)
    frames = 2:hist_len+1*framerate

    PolyPlanning.record(fig, "videos/dcol_$video_name.mp4", frames;
        framerate) do frame
        #@infiltrate
        if frame <= hist_len
            update_fig(dcol_sol.θ_history[frame])
        end
    end
end