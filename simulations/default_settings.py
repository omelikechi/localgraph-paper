# Default simulation settings for paper

default_settings = {}

default_settings['linear_sparse_n100'] = dict(
	n = 100,
	p = 200,
	snr = None,
	block_sizes = [1, 4, 195],
	block_degree = [0, 0, 2],
	connector_degree = [4, 6],
	block_magnitude = [1, 1, 1],
	connector_magnitude = [1, 1],
	lmin = 0.01,
	lmax = 10,
	radii = [1, 2, 3, 4],
	qpath_max = 0.2,
	fdr_local = [0.2, 0.1, 0.1, 0.1],
	fdr = 0.1,
	ipss_selector = 'l1',
	do_nonlinear = False
)

default_settings['nonlinear_sparse_n100'] = dict(
	n = 100,
	p = 200,
	snr = 4,
	block_sizes = [1, 4, 195],
	block_degree = [0, 0, 2],
	connector_degree = [4, 2],
	block_magnitude = [1, 1, 1],
	connector_magnitude = [1, 1],
	lmin = 0.01,
	lmax = 10,
	radii = [1, 2, 3, 4],
	qpath_max = 0.2,
	fdr_local = [0.2, 0.05, 0.05, 0.05],
	fdr = 0.1,
	ipss_selector = 'gb',
	do_nonlinear = True
)

default_settings['linear_dense_n100'] = dict(
	n = 100,
	p = 200,
	snr = None,
	block_sizes = [1, 20, 179],
	block_degree = [0, 5, 10],
	connector_degree = [20, 5],
	block_magnitude = [1, 1, 1],
	connector_magnitude = [1, 1],
	lmin = 0.01,
	lmax = 10,
	radii = [1, 2, 3],
	qpath_max = 0.5,
	fdr_local = [0.5, 0.15, 0.15],
	fdr = 0.1,
	ipss_selector = 'l1',
	do_nonlinear = False
)

default_settings['linear_dense_n500'] = dict(
	n = 500,
	p = 200,
	snr = None,
	block_sizes = [1, 20, 179],
	block_degree = [0, 5, 10],
	connector_degree = [20, 5],
	block_magnitude = [1, 1, 1],
	connector_magnitude = [1, 1],
	lmin = 0.01,
	lmax = 10,
	radii = [1, 2, 3],
	qpath_max = 0.2,
	fdr_local = [0.2, 0.1, 0.1, 0.1],
	fdr = 0.1,
	ipss_selector = 'l1',
	do_nonlinear = False
)

default_settings['nonlinear_dense_n500'] = dict(
	n = 500,
	p = 200,
	snr = 4,
	block_sizes = [1, 20, 179],
	block_degree = [0, 5, 10],
	connector_degree = [20, 5],
	block_magnitude = [1, 1, 1],
	connector_magnitude = [1, 1],
	lmin = 0.01,
	lmax = 10,
	radii = [1, 2, 3],
	qpath_max = 0.5,
	fdr_local = [0.2, 0.05, 0.05, 0.05],
	fdr = 0.1,
	ipss_selector = 'gb',
	do_nonlinear = True
)


