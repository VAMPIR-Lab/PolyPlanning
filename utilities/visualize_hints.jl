# visualize all
plot_n_rows = 5
plot_n_cols = 3
PolyPlanning.visualize_multi(x0s, maps, our_sols, T, ego_poly; n_rows=plot_n_rows, n_cols=plot_n_cols, type="nonsmooth")
PolyPlanning.visualize_multi(x0s, maps, sep_sols, T, ego_poly; n_rows=plot_n_rows, n_cols=plot_n_cols, type="sep_planes")
PolyPlanning.visualize_multi(x0s, maps, dcol_sols, T, ego_poly; n_rows=plot_n_rows, n_cols=plot_n_cols, type="dcol")
PolyPlanning.visualize_multi(x0s, maps, kkt_sols, T, ego_poly; n_rows=plot_n_rows, n_cols=plot_n_cols, type="direct_kkt")

# visualize successes
PolyPlanning.visualize_multi(x0s, maps, our_sols, our_bins.success, T, ego_poly; n_rows=plot_n_rows, n_cols=plot_n_cols, type="nonsmooth")
PolyPlanning.visualize_multi(x0s, maps, sep_sols, sep_bins.success, T, ego_poly; n_rows=plot_n_rows, n_cols=plot_n_cols, type="sep_planes")
PolyPlanning.visualize_multi(x0s, maps, dcol_sols, dcol_bins.success, T, ego_poly; n_rows=plot_n_rows, n_cols=plot_n_cols, type="dcol")
PolyPlanning.visualize_multi(x0s, maps, kkt_sols, kkt_bins.success, T, ego_poly; n_rows=plot_n_rows, n_cols=plot_n_cols, type="direct_kkt")

# visualize common success
our_sep_common_success = PolyPlanning.find_bin_common(our_bins.success, sep_bins.success)
PolyPlanning.visualize_multi(x0s, maps, our_sols, our_sep_common_success, T, ego_poly; n_rows=plot_n_rows, n_cols=plot_n_cols, type="nonsmooth")
PolyPlanning.visualize_multi(x0s, maps, kkt_sols, our_sep_common_success, T, ego_poly; n_rows=plot_n_rows, n_cols=plot_n_cols, type="sep_planes")

our_dcol_common_success = PolyPlanning.find_bin_common(our_bins.success, dcol_bins.success)
PolyPlanning.visualize_multi(x0s, maps, our_sols, our_dcol_common_success, T, ego_poly; n_rows=plot_n_rows, n_cols=plot_n_cols, type="nonsmooth")
PolyPlanning.visualize_multi(x0s, maps, dcol_sols, our_dcol_common_success, T, ego_poly; n_rows=plot_n_rows, n_cols=plot_n_cols, type="dcol")

our_kkt_common_success = PolyPlanning.find_bin_common(our_bins.success, kkt_bins.success)
PolyPlanning.visualize_multi(x0s, maps, our_sols, our_kkt_common_success, T, ego_poly; n_rows=plot_n_rows, n_cols=plot_n_cols, type="nonsmooth")
PolyPlanning.visualize_multi(x0s, maps, kkt_sols, our_kkt_common_success, T, ego_poly; n_rows=plot_n_rows, n_cols=plot_n_cols, type="direct_kkt")