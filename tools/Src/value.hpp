#include "headers.hpp"

template <typename T>
class value{
	public:
		int N;
		int ind;
		Eigen::Array<T, -1, 1> array;
		T val;
		
		value(Eigen::Array<T, -1, 1>, int);
		value(Eigen::Array<T, -1, 1>);
		value<T> operator+(value<T>);
		value<T> operator-(value<T>);
		value<T> operator-();
		
		template <typename Y>
		friend std::ostream& operator<<(std::ostream& os, value<Y>);
		void update_index(int);
};

template <typename T>
value<T>::value(Eigen::Array<T, -1, 1> arr, int i){
	N = arr.rows();
	array = arr;
	ind = i;
	val = array(ind);
}

template <typename T>
value<T>::value(Eigen::Array<T, -1, 1> arr){
	N = arr.rows();
	array = arr;
}

template <typename T>
value<T> value<T>::operator+(value<T> other){
	value<T> e(array, ind + other.ind - (N - 1)/2);
	return e;
}

template <typename T>
value<T> value<T>::operator-(value<T> other){
	
	return *this + (-other);
}

template <typename T>
value<T> value<T>::operator-(){
	//std::cout << "N: " << N << "ind: " << ind << "\n";
	//std::cout << "minus: " << N - 1 - ind << "\n";
	value<T> e(array, N - 1 - ind);
	return e;
}


template <typename T>
std::ostream& operator<<(std::ostream& os, value<T> other){
	os << other.val;
	return os;
}

template <typename T>
void value<T>::update_index(int i){
	/* If you want to change the value of 'ind' while updating the value 'val' use this function.
	 * This may not be needed because you may not need val at all. In that case, changing 'ind' by
	 * explicitly setting it to the new value is faster. e.g.: 
	 * value e(array,2); 
	 * e.ind = 3;
	 * */
	ind = i;
	val = array(i);
}
