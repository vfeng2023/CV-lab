// CPP program to demonstrate working of unordered_map
// for user defined data types.
#include <bits/stdc++.h>
using namespace std;

class Person {
	private:
	string first, last;
	public:
	Person(string f, string l)
	{
		first = f;
		last = l;
	}
	string getFirst(){
		return first;
	}
	string getLast(){
		return last;
	}
	bool operator==(const Person& p) const
	{
		return first == p.first && last == p.last;
	}
};

class MyHashFunction {
public:

	// We use predefined hash functions of strings
	// and define our hash function as XOR of the
	// hash values.
	size_t operator()(const Person& p) const
	{
		return (hash<string>()(p.first)) ^
			(hash<string>()(p.last));
	}
};

// Driver code
int main()
{
	unordered_map<Person, int, MyHashFunction> um;
	Person p1("kartik", "kapoor");
	Person p2("Ram", "Singh");
	Person p3("Laxman", "Prasad");

	um[p1] = 100;
	um[p2] = 200;
	um[p3] = 100;

	for (auto e : um) {
		cout << "[" << e.firstgetFirst() << ", "<< e.second.getLast()<< "] = > " << e.second << '\n';
	}

	return 0;
}
