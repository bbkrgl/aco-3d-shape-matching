template <class T, class V>
class FibonacciHeap;

template <class T, class V>
struct node {
private:
	node<T, V>* prev;
	node<T, V>* next;
	node<T, V>* child;
	node<T, V>* parent;

	T key;
	V value;
	int degree;
	bool marked;

public:
	friend class FibonacciHeap<T, V>;

	node<T, V>* getPrev() { return prev; }
	node<T, V>* getNext() { return next; }
	node<T, V>* getChild() { return child; }
	node<T, V>* getParent() { return parent; }

	V getValue() { return value; }
	bool isMarked() { return marked; }

	bool hasChildren() { return child; }
	bool hasParent() { return parent; }
};

template <class T, class V>
class FibonacciHeap {
protected:
	node<T, V>* heap;

public:
	FibonacciHeap() { heap = _empty(); }

	virtual ~FibonacciHeap()
	{
		if (heap)
			_deleteAll(heap);
	}

	node<T, V>* insert(T key, V value)
	{
		node<T, V>* ret = _singleton(key, value);
		heap = _merge(heap, ret);
		return ret;
	}

	void merge(FibonacciHeap& other)
	{
		heap = _merge(heap, other.heap);
		other.heap = _empty();
	}

	bool isEmpty() { return heap == nullptr; }

	T getMinimum() { return heap->key; }

	V removeMinimum()
	{
		node<T, V>* old = heap;
		heap = _removeMinimum(heap);
		V ret = old->value;
		delete old;
		return ret;
	}

	void decreaseKey(node<T, V>* n, V value) { heap = _decreaseKey(heap, n, value); }

	node<T, V>* find(T key) { return _find(heap, key); }

private:
	node<T, V>* _empty() { return nullptr; }

	node<T, V>* _singleton(T key, V value)
	{
		node<T, V>* n = new node<T, V>;
		n->key = key;
		n->value = value;
		n->prev = n->next = n;
		n->degree = 0;
		n->marked = false;
		n->child = nullptr;
		n->parent = nullptr;
		return n;
	}

	node<T, V>* _merge(node<T, V>* a, node<T, V>* b)
	{
		if (a == nullptr)
			return b;

		if (b == nullptr)
			return a;

		if (a->value > b->value) {
			node<T, V>* temp = a;
			a = b;
			b = temp;
		}

		node<T, V>* an = a->next;
		node<T, V>* bp = b->prev;
		a->next = b;
		b->prev = a;
		an->prev = bp;
		bp->next = an;
		return a;
	}

	void _deleteAll(node<T, V>* n)
	{
		if (n != nullptr) {
			node<T, V>* c = n;
			do {
				node<T, V>* d = c;
				c = c->next;
				_deleteAll(d->child);
				delete d;
			} while (c !=  n);
		}
	}

	void _addChild(node<T, V>* parent, node<T, V>* child)
	{
		child->prev = child->next = child;
		child->parent = parent;
		parent->degree++;
		parent->child = _merge(parent->child, child);
	}

	void _unMarkAndUnParentAll(node<T, V>* n)
	{
		if (n == nullptr)
			return;

		node<T, V>* c = n;
		do {
			c->marked = false;
			c->parent = nullptr;
			c = c->next;
		} while (c !=  n);
	}

	node<T, V>* _removeMinimum(node<T, V>* n)
	{
		_unMarkAndUnParentAll(n->child);

		if (n->next == n) {
			n = n->child;
		} else {
			n->next->prev = n->prev;
			n->prev->next = n->next;
			n = _merge(n->next, n->child);
		}

		if (n == nullptr)
			return n;

		node<T, V>* trees[64] = {nullptr};
		while (true) {
			if (trees[n->degree] !=  nullptr) {
				node<T, V>* t = trees[n->degree];

				if (t == n)
					break;

				trees[n->degree] = nullptr;

				if (n->value < t->value) {
					t->prev->next = t->next;
					t->next->prev = t->prev;
					_addChild(n, t);
				} else {
					t->prev->next = t->next;
					t->next->prev = t->prev;
					if (n->next == n) {
						t->next = t->prev = t;
						_addChild(t, n);
						n = t;
					} else {
						n->prev->next = t;
						n->next->prev = t;
						t->next = n->next;
						t->prev = n->prev;
						_addChild(t, n);
						n = t;
					}
				}

				continue;
			} else {
				trees[n->degree] = n;
			}

			n = n->next;
		}

		node<T, V>* min = n;
		node<T, V>* start = n;
		do {
			if (n->value < min->value)
				min = n;

			n = n->next;
		} while (n !=  start);

		return min;
	}

	node<T, V>* _cut(node<T, V>* heap, node<T, V>* n)
	{
		if (n->next == n) {
			n->parent->child = nullptr;
		} else {
			n->next->prev = n->prev;
			n->prev->next = n->next;
			n->parent->child = n->next;
		}

		n->next = n->prev = n;
		n->marked = false;
		return _merge(heap, n);
	}

	node<T, V>* _decreaseKey(node<T, V>* heap, node<T, V>* n, V value)
	{
		if (n->value < value)
			return heap;

		n->value = value;
		if (n->parent && n->value < n->parent->value) {
			heap = _cut(heap, n);
			node<T, V>* parent = n->parent;
			n->parent = nullptr;
			while (parent !=  nullptr && parent->marked) {
				heap = _cut(heap, parent);
				n = parent;
				parent = n->parent;
				n->parent = nullptr;
			}

			if (parent !=  nullptr && parent->parent !=  nullptr)
				parent->marked = true;
		} else if (n->value < heap->value) {
			heap = n;
		}

		return heap;
	}

	node<T, V>* _find(node<T, V>* heap, T key)
	{
		node<T, V>* n = heap;

		if (n == nullptr)
			return nullptr;

		do {
			if (n->key == key)
				return n;

			node<T, V>* ret = _find(n->child, key);
			if (ret)
				return ret;

			n = n->next;
		} while (n !=  heap);

		return nullptr;
	}
};