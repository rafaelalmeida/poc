#ifndef COUNTER_H
#define COUNTER_H

// Helper class to count things
template <typename T>
class Counter {
	std::map<T, int> _map;

	public:
		void inc(T item);
		T top();
		std::map<T, int> getCounts();
		int getCount(T item);
};

template <typename T>
void Counter<T>::inc(T item) {
	if (_map.count(item) == 0) {
		_map[item] = 0;
	}

	_map[item]++;
}

template <typename T> 
T Counter<T>::top() {
	int maxVal = 0;
	T cMax;

	for (auto i : _map) {
		if (i.second > maxVal) {
			maxVal = i.second;
			cMax = i.first;
		}
	}

	return cMax;
}

template <typename T>
std::map<T, int> Counter<T>::getCounts() {
	return _map;
}

template <typename T>
int Counter<T>::getCount(T item) {
	auto val = _map.find(item);
	if (val == _map.end()) {
		return 0;
	}
	else {
		return val->second;
	}
}

#endif // COUNTER_H
