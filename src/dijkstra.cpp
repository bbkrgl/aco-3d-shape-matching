#include "dijkstra.h"

void get_path(Eigen::MatrixXi prev, std::vector<std::vector<int>> &edges, int d)
{
	while(prev(d) >= 0) {
		std::vector<int> v;
		v.push_back(d);
		v.push_back(prev(d));
		edges.push_back(v);
		d = prev(d);
	}
}

void dijkstra(int p, std::vector<std::vector<int>> &adj_list, Eigen::MatrixXd &V, Eigen::VectorXd &Q, Eigen::VectorXi &prev)
{
	prev.fill(-1);
	Q.fill(DBL_MAX);
	Q(p) = 0;

	MHeap mh(adj_list.size());
	mh.insertKey(p, 0);

	for (int i = 0; i < adj_list.size(); i++) {
		if (i == p)
			continue;

		mh.insertKey(i, DBL_MAX);
	}

	while(!mh.isEmpty()) {
		int u = mh.getMin();
		mh.extractMin();
		for (int i = 0; i < adj_list[u].size(); i++) {
			int s = u;
			int d = adj_list[u][i];

			double distance = 
				sqrt(pow(V(s,0) - V(d,0), 2) + pow(V(s,1) - V(d,1), 2) + pow(V(s,2) - V(d,2), 2));

			if (Q(d) > distance + Q(s)) {
				Q(d) = distance + Q(s);
				prev(d) = s;
				mh.decreaseKey(d, distance + Q(s));
			}
		}
	}
}