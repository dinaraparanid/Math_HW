#include <iostream>
#include <iomanip>
#include <vector>
#include <optional>
#include <variant>
#include <functional>
#include <iterator>
#include <cmath>
#include <numeric>

// #define DEBUG

namespace agla {
	namespace mtx {
		template <typename T> struct matrix {

			// ########################## Matrix Row ##########################

			class matrix_row {
				std::vector<T> row;

			public:

				// ----------------------- Iterators -----------------------

				class const_iterator;

				// ########################## Iterator ##########################

				class iterator {
				public:
					using iterator_category = std::random_access_iterator_tag;
					using difference_type = std::ptrdiff_t;
					using value_type = T;
					using pointer = value_type*;
					using reference = value_type&;

				private:
					friend class matrix_row;
					std::vector<T>::iterator row_it;
					explicit iterator(std::vector<T>::iterator row_it) noexcept : row_it(row_it) {}

				public:
					~iterator() noexcept = default;

					// --------------- Dereference operators ---------------

					inline reference operator*() const noexcept { return *row_it; }
					inline pointer operator->() const noexcept { return &*row_it; }

					// --------------- Comparison operators ---------------

					[[nodiscard]] inline bool operator==(const iterator& other) const noexcept { return row_it == other.row_it; };
					[[nodiscard]] inline bool operator!=(const iterator& other) const noexcept { return row_it != other.row_it; };

					[[nodiscard]] inline bool operator==(const const_iterator& other) const noexcept { return row_it == other.row_it; };
					[[nodiscard]] inline bool operator!=(const const_iterator& other) const noexcept { return row_it != other.row_it; };

					// --------------- Movement operators ---------------

					inline iterator operator++() noexcept { ++row_it; return *this; }
					inline iterator operator--() noexcept { --row_it; return *this; }

					[[nodiscard]] inline iterator operator+(const difference_type move) const noexcept { return iterator(row_it + move); }
					[[nodiscard]] inline iterator operator-(const difference_type move) const noexcept { return iterator(row_it - move); }

					[[nodiscard]] inline difference_type operator-(const const_iterator it) const noexcept { return row_it - it.row_it; }
					[[nodiscard]] inline difference_type operator-(const iterator it) const noexcept { return row_it - it.row_it; }

					[[nodiscard]] inline iterator operator+=(const difference_type move) noexcept {
						row_it += move;
						return *this;
					}

					[[nodiscard]] inline iterator operator-=(const difference_type move) noexcept {
						row_it -= move;
						return *this;
					}
				};

				// ########################## Const Iterator ##########################

				class const_iterator {
				public:
					using iterator_category = std::random_access_iterator_tag;
					using difference_type = std::ptrdiff_t;
					using value_type = T;
					using pointer = const value_type*;
					using reference = const value_type&;

				private:
					friend class matrix_row;
					std::vector<T>::const_iterator row_it;
					explicit const_iterator(std::vector<T>::const_iterator row_it) noexcept : row_it(row_it) {}

				public:
					~const_iterator() noexcept = default;

					// --------------- Dereference operators ---------------

					inline reference operator*() const noexcept { return *row_it; }
					inline pointer operator->() const noexcept { return &*row_it; }

					// --------------- Comparison operators ---------------

					[[nodiscard]] inline bool operator==(const iterator& other) const noexcept { return row_it == other.row_it; };
					[[nodiscard]] inline bool operator!=(const iterator& other) const noexcept { return row_it != other.row_it; };

					[[nodiscard]] inline bool operator==(const const_iterator& other) const noexcept { return row_it == other.row_it; };
					[[nodiscard]] inline bool operator!=(const const_iterator& other) const noexcept { return row_it != other.row_it; };

					// --------------- Movement operators ---------------

					inline const_iterator operator++() noexcept { ++row_it; return *this; }
					inline const_iterator operator--() noexcept { --row_it; return *this; }

					[[nodiscard]] inline const_iterator operator+(const difference_type move) const noexcept { return const_iterator(row_it + move); }
					[[nodiscard]] inline const_iterator operator-(const difference_type move) const noexcept { return const_iterator(row_it - move); }

					[[nodiscard]] inline difference_type operator-(const const_iterator it) const noexcept { return row_it - it.row_it; }
					[[nodiscard]] inline difference_type operator-(const iterator it) const noexcept { return row_it - it.row_it; }

					[[nodiscard]] inline const_iterator operator+=(const difference_type move) noexcept {
						row_it += move;
						return *this;
					}

					[[nodiscard]] inline const_iterator operator-=(const difference_type move) noexcept {
						row_it -= move;
						return *this;
					}
				};

				// ----------------------- Constructors -----------------------

				matrix_row() noexcept = default;

				explicit matrix_row(const std::size_t size) noexcept :
						row(std::vector<T>(size)) {}

				matrix_row(const std::size_t size, const T elem) noexcept :
						row(std::vector<T>(size, elem)) {}

				explicit matrix_row(const std::vector<T>& row) noexcept : row(row) {}
				explicit matrix_row(std::vector<T>&& row) noexcept : row(row) {}

				// ----------------------- Accessors -----------------------

				[[nodiscard]] inline std::size_t size() const noexcept { return row.size(); }

				[[nodiscard]] inline T& get_unchecked(const std::size_t index) noexcept {
					return row[index];
				}

				[[nodiscard]] inline const T& get_unchecked(const std::size_t index) const noexcept {
					return row[index];
				}

				[[nodiscard]] inline std::optional<std::reference_wrapper<T>> operator[](const std::size_t index) noexcept {
					if (index >= row.size()) return std::nullopt;
					return std::make_optional(std::ref(get_unchecked(index)));
				}

				[[nodiscard]] inline std::optional<matrix_row> operator+(const matrix_row& other) const noexcept {
					const auto sz = size();

					if (sz != other.size())
						return std::nullopt;

					matrix_row result(sz);

					for (std::size_t i = 0; i < sz; ++i)
						result.get_unchecked(i) = get_unchecked(i) + other.get_unchecked(i);

					return std::make_optional(result);
				}

				[[nodiscard]] inline std::optional<matrix_row> operator-(const matrix_row& other) const noexcept {
					const auto sz = size();

					if (sz != other.size())
						return std::nullopt;

					matrix_row result(sz);

					for (std::size_t i = 0; i < sz; ++i)
						result.get_unchecked(i) = get_unchecked(i) - other.get_unchecked(i);

					return std::make_optional(result);
				}

				// ----------------------- Iterators -----------------------

				[[nodiscard]] inline iterator begin() noexcept { return iterator(row.begin()); }
				[[nodiscard]] inline const_iterator begin() const noexcept { return const_iterator(row.begin()); }

				[[nodiscard]] inline iterator end() noexcept { return iterator(row.end()); }
				[[nodiscard]] inline const_iterator end() const noexcept { return const_iterator(row.end()); }
			};

		protected:
			std::vector<matrix_row> mtx;

			// ----------------------- Row Iterators -----------------------

			class const_row_iterator;

			// ########################## Row Iterator ##########################

			class row_iterator {
			public:
				using iterator_category = std::random_access_iterator_tag;
				using difference_type = std::ptrdiff_t;
				using value_type = matrix_row;
				using pointer = value_type*;
				using reference = value_type&;

			private:
				friend class matrix;
				std::vector<matrix_row>::iterator it;
				explicit row_iterator(std::vector<matrix_row>::iterator it) noexcept : it(it) {}

			public:
				~row_iterator() noexcept = default;

				// --------------- Dereference operators ---------------

				inline reference operator*() const noexcept { return *it; }
				inline pointer operator->() const noexcept { return &*it; }

				// --------------- Comparison operators ---------------

				[[nodiscard]] inline bool operator==(const row_iterator& other) const noexcept { return it == other.it; };
				[[nodiscard]] inline bool operator!=(const row_iterator& other) const noexcept { return it != other.it; };

				[[nodiscard]] inline bool operator==(const const_row_iterator& other) const noexcept { return it == other.it; };
				[[nodiscard]] inline bool operator!=(const const_row_iterator& other) const noexcept { return it != other.it; };

				// --------------- Movement operators ---------------

				inline row_iterator operator++() noexcept { ++it; return *this; }
				inline row_iterator operator--() noexcept { --it; return *this; }

				[[nodiscard]] inline row_iterator operator+(const difference_type move) const noexcept { return row_iterator(it + move); }
				[[nodiscard]] inline row_iterator operator-(const difference_type move) const noexcept { return row_iterator(it - move); }

				[[nodiscard]] inline difference_type operator-(const row_iterator iter) const noexcept { return it - iter.it; }
				[[nodiscard]] inline difference_type operator-(const const_row_iterator iter) const noexcept { return it - iter.it; }

				[[nodiscard]] inline row_iterator operator+=(const difference_type move) noexcept {
					it += move;
					return *this;
				}

				[[nodiscard]] inline row_iterator operator-=(const difference_type move) noexcept {
					it -= move;
					return *this;
				}
			};

			// ########################## Const Iterator ##########################

			class const_row_iterator {
			public:
				using iterator_category = std::random_access_iterator_tag;
				using difference_type = std::ptrdiff_t;
				using value_type = matrix_row;
				using pointer = const value_type*;
				using reference = const value_type&;

			private:
				friend class matrix;
				std::vector<matrix_row>::const_iterator it;
				explicit const_row_iterator(std::vector<matrix_row>::const_iterator it) noexcept : it(it) {}

			public:
				~const_row_iterator() noexcept = default;

				// --------------- Dereference operators ---------------

				inline reference operator*() const noexcept { return *it; }
				inline pointer operator->() const noexcept { return &*it; }

				// --------------- Comparison operators ---------------

				[[nodiscard]] inline bool operator==(const row_iterator& other) const noexcept { return it == other.it; };
				[[nodiscard]] inline bool operator!=(const row_iterator& other) const noexcept { return it != other.it; };

				[[nodiscard]] inline bool operator==(const const_row_iterator& other) const noexcept { return it == other.it; };
				[[nodiscard]] inline bool operator!=(const const_row_iterator& other) const noexcept { return it != other.it; };

				// --------------- Movement operators ---------------

				inline const_row_iterator operator++() noexcept { ++it; return *this; }
				inline const_row_iterator operator--() noexcept { --it; return *this; }

				[[nodiscard]] inline const_row_iterator operator+(const difference_type move) const noexcept { return const_row_iterator(it + move); }
				[[nodiscard]] inline const_row_iterator operator-(const difference_type move) const noexcept { return const_row_iterator(it - move); }

				[[nodiscard]] inline difference_type operator-(const row_iterator iter) const noexcept { return it - iter.it; }
				[[nodiscard]] inline difference_type operator-(const const_row_iterator iter) const noexcept { return it - iter.it; }

				[[nodiscard]] inline const_row_iterator operator+=(const difference_type move) noexcept {
					it += move;
					return *this;
				}

				[[nodiscard]] inline const_row_iterator operator-=(const difference_type move) noexcept {
					it -= move;
					return *this;
				}
			};

		public:

			// ----------------------- Iterators -----------------------

			class const_iterator;

			// ########################## Iterator ##########################

			class iterator {
			public:
				using iterator_category = std::bidirectional_iterator_tag;
				using difference_type = std::ptrdiff_t;
				using value_type = T;
				using pointer = value_type*;
				using reference = value_type&;

			private:
				friend class matrix;
				row_iterator row_it;
				matrix_row::iterator it;

				const row_iterator rows_begin;
				const row_iterator rows_end;

				iterator(
						row_iterator row_it,
						matrix_row::iterator it,
						const row_iterator rows_begin,
						const row_iterator rows_end
				) noexcept : row_it(row_it), it(it), rows_begin(rows_begin), rows_end(rows_end) {}

			public:
				~iterator() noexcept = default;

				// --------------- Dereference operators ---------------

				inline reference operator*() const noexcept { return *it; }
				inline pointer operator->() const noexcept { return &*it; }

				// --------------- Comparison operators ---------------

				[[nodiscard]] inline bool operator==(const iterator& other) const noexcept { return it == other.it; };
				[[nodiscard]] inline bool operator!=(const iterator& other) const noexcept { return it != other.it; };

				[[nodiscard]] inline bool operator==(const const_iterator& other) const noexcept { return it == other.it; };
				[[nodiscard]] inline bool operator!=(const const_iterator& other) const noexcept { return it != other.it; };

				// --------------- Movement operators ---------------

				inline iterator operator++() noexcept {
					const auto row_end_it = row_it->end();

					if (std::next(it) == row_end_it)
						it = (row_it + 1 == rows_end) ? row_end_it : (++row_it)->begin();
					else
						++it;

					return *this;
				}

				inline iterator operator--() noexcept {
					const auto row_begin_it = row_it->begin();

					if (it == row_begin_it) {
						if (row_it != rows_begin)
							it = std::prev((--row_it)->end());
					} else { --it; }

					return *this;
				}
			};

			[[nodiscard]] inline iterator iter(row_iterator row_it, matrix_row::iterator it) noexcept {
				return iterator(row_it, it, rows_begin(), rows_end());
			}

			// ########################## Const Iterator ##########################

			class const_iterator {
			public:
				using iterator_category = std::bidirectional_iterator_tag;
				using difference_type = std::ptrdiff_t;
				using value_type = T;
				using pointer = const value_type*;
				using reference = const value_type&;

			private:
				friend class matrix;
				const_row_iterator row_it;
				matrix_row::const_iterator it;

				const const_row_iterator rows_begin;
				const const_row_iterator rows_end;

				const_iterator(
						const_row_iterator row_it,
						matrix_row::const_iterator it,
						const const_row_iterator rows_begin,
						const const_row_iterator rows_end
				) noexcept : row_it(row_it), it(it), rows_begin(rows_begin), rows_end(rows_end) {}

			public:
				~const_iterator() noexcept = default;

				// --------------- Dereference operators ---------------

				inline reference operator*() const noexcept { return *it; }
				inline pointer operator->() const noexcept { return &*it; }

				// --------------- Comparison operators ---------------

				[[nodiscard]] inline bool operator==(const iterator& other) const noexcept { return it == other.it; };
				[[nodiscard]] inline bool operator!=(const iterator& other) const noexcept { return it != other.it; };

				[[nodiscard]] inline bool operator==(const const_iterator& other) const noexcept { return it == other.it; };
				[[nodiscard]] inline bool operator!=(const const_iterator& other) const noexcept { return it != other.it; };

				// --------------- Movement operators ---------------

				inline const_iterator operator++() noexcept {
					const auto row_end_it = row_it->end();

					if (std::next(it) == row_end_it)
						it = (row_it + 1 == rows_end) ? row_end_it : (++row_it)->begin();
					else
						++it;

					return *this;
				}

				inline const_iterator operator--() noexcept {
					const auto row_begin_it = row_it->begin();

					if (it == row_begin_it) {
						if (row_it != rows_begin)
							it = std::prev((--row_it)->end());
					} else { --it; }

					return *this;
				}
			};

			[[nodiscard]] inline const_iterator const_iter(
					const_row_iterator row_it,
					matrix_row::const_iterator it
			) const noexcept {
				return const_iterator(row_it, it, rows_begin(), rows_end());
			}

			// ----------------------- Constructors -----------------------

			explicit matrix(const std::size_t size) noexcept :
					mtx(std::vector<matrix_row>(size, matrix_row(size))) {}

			matrix(const std::size_t rows, const std::size_t columns) noexcept :
					mtx(std::vector<matrix_row>(rows, matrix_row(columns))) {}

			matrix(const std::size_t rows, const std::vector<T>& row) noexcept :
					mtx(std::vector<matrix_row>(rows, matrix_row(row))) {}

			explicit matrix(const std::vector<std::vector<T>>& matrix) noexcept {
				for (const auto& row : matrix)
					mtx.emplace_back(row);
			}

			explicit matrix(std::vector<std::vector<T>>&& matrix) noexcept {
				for (auto&& row : matrix)
					mtx.emplace_back(row);
			}

			~matrix() noexcept = default;

			// ----------------------- Accessors -----------------------

			[[nodiscard]] inline std::size_t rows_number() const noexcept { return mtx.size(); }
			[[nodiscard]] inline std::size_t columns_number() const noexcept { return mtx.front().size(); }

			[[nodiscard]] inline matrix_row& get_unchecked(const std::size_t index) noexcept {
				return mtx[index];
			}

			[[nodiscard]] const matrix_row& get_unchecked(const std::size_t index) const noexcept {
				return mtx[index];
			}

			[[nodiscard]] virtual inline std::optional<std::reference_wrapper<matrix_row>> operator[](const std::size_t index) noexcept {
				if (index >= mtx.size()) return std::nullopt;
				return std::make_optional(std::ref(get_unchecked(index)));
			}

			// ----------------------- Operations -----------------------

			[[nodiscard]] inline std::optional<matrix> operator+(const matrix& other) const noexcept {
				const auto rows = rows_number();
				const auto columns = columns_number();

				if (rows != other.rows_number() || columns != other.columns_number())
					return std::nullopt;

				matrix result(rows, columns);

				auto this_it = begin();
				auto other_it = other.begin();

				for (auto& elem : result) {
					elem = *this_it + *other_it;
					++this_it, ++other_it;
				}

				return std::make_optional(result);
			}

			[[nodiscard]] inline std::optional<matrix> operator-(const matrix& other) const noexcept {
				const auto rows = rows_number();
				const auto columns = columns_number();

				if (rows != other.rows_number() || columns != other.columns_number())
					return std::nullopt;

				matrix result(rows, columns);

				auto this_it = begin();
				auto other_it = other.begin();

				for (auto& elem : result) {
					elem = *this_it - *other_it;
					++this_it, ++other_it;
				}

				return std::make_optional(result);
			}

			[[nodiscard]] inline std::optional<matrix> operator* (const matrix& other) const noexcept {
				if (columns_number() != other.rows_number())
					return std::nullopt;

				matrix result(rows_number(), other.columns_number());

				for (int i = 0; i < result.rows_number(); ++i)
					for (int q = 0; q < result.columns_number(); ++q)
						for (int j = 0; j < columns_number(); ++j)
							result.get_unchecked(i).get_unchecked(q) +=
									this->get_unchecked(i).get_unchecked(j) *
											other.get_unchecked(j).get_unchecked(q);

				return std::make_optional(result);
			}

			[[nodiscard]] inline bool operator== (const matrix& other) const noexcept {
				if (rows_number() != other.rows_number() || columns_number() != other.columns_number())
					return false;

				return std::equal(begin(), end(), other.begin());
			}

			[[nodiscard]] inline bool operator!= (const matrix& other) const noexcept {
				return !(*this == other);
			}

			[[nodiscard]] inline matrix transposed() const noexcept {
				matrix result(columns_number(), rows_number());

				for (int i = 0; i < result.rows_number(); ++i)
					for (int q = 0; q <result.columns_number(); ++q)
						result.get_unchecked(i).get_unchecked(q) = get_unchecked(q).get_unchecked(i);

				return result;
			}

			inline matrix& operator=(const matrix& matrix) noexcept {
				mtx = matrix.mtx;
				return *this;
			}

			#ifdef DEBUG
			inline void debug_mtx(std::size_t& step, const char* const msg) const noexcept {
				std::printf("step %zu: %s\n", step, msg);
				std::cout << *this;
			}
			#endif

			// ----------------------- Iterators -----------------------

			[[nodiscard]] inline row_iterator rows_begin() noexcept { return row_iterator(mtx.begin()); }
			[[nodiscard]] inline const_row_iterator rows_begin() const noexcept { return const_row_iterator(mtx.begin()); }

			[[nodiscard]] inline row_iterator rows_end() noexcept { return row_iterator(mtx.end()); }
			[[nodiscard]] inline const_row_iterator rows_end() const noexcept { return const_row_iterator(mtx.end()); }

			[[nodiscard]] inline iterator begin() noexcept { return iter(rows_begin(), mtx.front().begin()); }
			[[nodiscard]] inline const_iterator begin() const noexcept { return const_iter(rows_begin(), mtx.front().begin()); }

			[[nodiscard]] inline iterator end() noexcept { return iter(rows_end(), mtx.back().end()); }
			[[nodiscard]] inline const_iterator end() const noexcept { return const_iter(rows_end(), mtx.back().end()); }
		};

		template <typename T> class square_matrix : public matrix<T> {
			explicit square_matrix(const matrix<T>& mtx) noexcept : matrix<T>(mtx) {}

		public:
			explicit square_matrix(const std::size_t size) noexcept : matrix<T>(size) {}

			square_matrix(const std::size_t size, const T& elem) noexcept : matrix<T>(size, std::vector<T>(size, elem)) {}
			square_matrix(const std::size_t size, T&& elem) noexcept : matrix<T>(size, std::vector<T>(size, elem)) {}

			square_matrix(const std::size_t rows, const std::vector<T>& row) noexcept : matrix<T>(rows, row) {
				if (rows != row.size())
					throw std::length_error("Number of rows does not equal to number of columns");
			}

			explicit square_matrix(const std::vector<std::vector<T>>& mtx) noexcept : matrix<T>(mtx) {}
			explicit square_matrix(std::vector<std::vector<T>>&& mtx) noexcept : matrix<T>(mtx) {}

			~square_matrix() noexcept = default;

			[[nodiscard]] inline std::size_t size() const noexcept {
				return this->rows_number();
			}

			[[nodiscard]] inline std::optional<square_matrix> operator+(const square_matrix& other) const noexcept {
				auto res = *static_cast<const matrix<T>*>(this) + static_cast<matrix<T>>(other);
				return res.has_value() ? std::make_optional(static_cast<square_matrix<T>>(res.value())) : std::nullopt;
			}

			[[nodiscard]] inline std::optional<square_matrix> operator-(const square_matrix& other) const noexcept {
				auto res = *static_cast<const matrix<T>*>(this) - static_cast<matrix<T>>(other);
				return res.has_value() ? std::make_optional(static_cast<square_matrix<T>>(res.value())) : std::nullopt;
			}

			inline square_matrix& operator=(const square_matrix& matrix) noexcept {
				this->mtx = matrix.mtx;
				return *this;
			}

			explicit square_matrix(matrix<T>&& mtx) noexcept : matrix<T>(mtx) {}
		public:
			[[nodiscard]] inline T determinant() const noexcept {
				auto copy = *this;
				const auto size = copy.size();
				T acc = 1;

				#ifdef DEBUG
				std::size_t step = 1;
				#endif

				for (std::size_t i = 0; i < size; ++i) {
					auto diag = copy.get_unchecked(i).get_unchecked(i);
					auto diag_index = i;

					for (std::size_t q = i + 1; q < size; ++q) {
						const auto cur_diag = copy.get_unchecked(q).get_unchecked(i);

						if (std::abs(cur_diag) > std::abs(diag)) {
							diag = cur_diag;
							diag_index = q;
						}
					}

					if (diag_index != i) {
						std::swap(copy.get_unchecked(i), copy.get_unchecked(diag_index));
						acc = -acc;

						#ifdef DEBUG
						copy.debug_mtx(step, "permutation");
						#endif
					}

					if (diag == 0) {
						#ifdef DEBUG
						std::puts("result:");
						#endif
						return 0;
					}

					for (std::size_t q = i + 1; q < size; ++q) {
						const auto ratio = copy.get_unchecked(q).get_unchecked(i) / diag;
						if (ratio == 0) continue;

						for (std::size_t k = 0; k < size; ++k)
							copy.get_unchecked(q).get_unchecked(k) -=
									copy.get_unchecked(i).get_unchecked(k) * ratio;

						#ifdef DEBUG
						copy.debug_mtx(step, "elimination");
						#endif
					}
				}

				for (std::size_t i = 0; i < size; ++i)
					acc *= copy.get_unchecked(i).get_unchecked(i);

				#ifdef DEBUG
				std::puts("result:");
				#endif

				return acc;
			}

			[[nodiscard]] inline std::optional<square_matrix<T>> inversed() const noexcept {
				const auto size = this->size();
				const auto aug_columns_num = 2 * size;
				matrix<T> aug_mtx(size, std::vector<T>(aug_columns_num, 0));

				for (std::size_t i = 0; i < size; ++i) {
					for (std::size_t q = 0; q < size; ++q)
						aug_mtx.get_unchecked(i).get_unchecked(q) = this->get_unchecked(i).get_unchecked(q);
					aug_mtx.get_unchecked(i).get_unchecked(i + size) = 1;
				}

				#ifdef DEBUG
				std::size_t step = 0;
				aug_mtx.debug_mtx(step, "Augmented Matrix");
				std::puts("Direct way:");
				#endif

				for (std::size_t i = 0; i < size; ++i) {
					auto diag = aug_mtx.get_unchecked(i).get_unchecked(i);
					auto diag_index = i;

					for (std::size_t q = i + 1; q < size; ++q) {
						const auto cur_diag = aug_mtx.get_unchecked(q).get_unchecked(i);

						if (std::abs(cur_diag) > std::abs(diag)) {
							diag = cur_diag;
							diag_index = q;
						}
					}

					if (diag_index != i) {
						std::swap(aug_mtx.get_unchecked(i), aug_mtx.get_unchecked(diag_index));
						#ifdef DEBUG
						aug_mtx.debug_mtx(step, "permutation");
						#endif
					}

					if (diag == 0) {
						#ifdef DEBUG
						std::puts("result:");
						#endif
						return std::nullopt;
					}

					for (std::size_t q = i + 1; q < size; ++q) {
						const auto ratio = aug_mtx.get_unchecked(q).get_unchecked(i) / diag;
						if (ratio == 0) continue;

						for (std::size_t k = 0; k < aug_columns_num; ++k)
							aug_mtx.get_unchecked(q).get_unchecked(k) -=
									aug_mtx.get_unchecked(i).get_unchecked(k) * ratio;

						#ifdef DEBUG
						aug_mtx.debug_mtx(step, "elimination");
						#endif
					}
				}

				#ifdef DEBUG
				std::puts("Way back:");
				#endif

				for (int i = size - 1; i >= 0; --i) {
					const auto diag = aug_mtx.get_unchecked(i).get_unchecked(i);

					for (int q = i - 1; q >= 0; --q) {
						const auto ratio = aug_mtx.get_unchecked(q).get_unchecked(i) / diag;
						if (ratio == 0) continue;

						for (std::size_t k = i; k < aug_columns_num; ++k)
							aug_mtx.get_unchecked(q).get_unchecked(k) -=
									aug_mtx.get_unchecked(i).get_unchecked(k) * ratio;

						#ifdef DEBUG
						aug_mtx.debug_mtx(step, "elimination");
						#endif
					}
				}

				for (std::size_t i = 0; i < size; ++i) {
					for (std::size_t q = size; q < aug_columns_num; ++q)
						aug_mtx.get_unchecked(i).get_unchecked(q) /=
								aug_mtx.get_unchecked(i).get_unchecked(i);
					aug_mtx.get_unchecked(i).get_unchecked(i) = 1;
				}

				#ifdef DEBUG
				std::puts("Diagonal normalization:");
				std::cout << aug_mtx;
				#endif

				square_matrix<T> result(size);

				for (std::size_t i = 0; i < size; ++i)
					for (std::size_t q = 0, qg = size; q < size; ++q, ++qg)
						result.get_unchecked(i).get_unchecked(q) =
								aug_mtx.get_unchecked(i).get_unchecked(qg);

				#ifdef DEBUG
				std::puts("result:");
				#endif

				return { std::move(result) };
			}
		};

		template <typename T> class identity_matrix : public square_matrix<T> {
			explicit identity_matrix(const square_matrix<T>& mtx) noexcept : square_matrix<T>(mtx) {}
			explicit identity_matrix(square_matrix<T>&& mtx) noexcept : square_matrix<T>(mtx) {}

		public:
			explicit identity_matrix(const std::size_t size) noexcept : square_matrix<T>(size, std::vector<T>(size, 0)) {
				for (std::size_t i = 0; i < size; ++i)
					this->get_unchecked(i).get_unchecked(i) = 1;
			}

			~identity_matrix() noexcept = default;

			inline identity_matrix& operator=(const identity_matrix& matrix) noexcept {
				this->mtx = matrix.mtx;
				return *this;
			}

		protected:
			[[nodiscard]] inline matrix<T>::matrix_row& get_unchecked(const std::size_t index) noexcept {
				return this->mtx[index];
			}
		};

		template <typename T> class elimination_matrix : public identity_matrix<T> {
			explicit elimination_matrix(const identity_matrix<T>& mtx) noexcept : identity_matrix<T>(mtx) {}
			explicit elimination_matrix(identity_matrix<T>&& mtx) noexcept : identity_matrix<T>(mtx) {}

		public:
			explicit elimination_matrix(
					const square_matrix<T>& mtx,
					const std::size_t row_ind,
					const std::size_t column_ind
			) noexcept : identity_matrix<T>(mtx.rows_number()) {
				const auto diag = mtx.get_unchecked(column_ind).get_unchecked(column_ind);
				const auto elem = mtx.get_unchecked(row_ind).get_unchecked(column_ind);
				this->get_unchecked(row_ind).get_unchecked(column_ind) = -elem / diag;
			}

			~elimination_matrix() noexcept = default;

			inline elimination_matrix& operator=(const elimination_matrix& matrix) noexcept {
				this->mtx = matrix.mtx;
				return *this;
			}
		};

		template <typename T> class permutation_matrix : public identity_matrix<T> {
			explicit permutation_matrix(const identity_matrix<T>& mtx) noexcept : identity_matrix<T>(mtx) {}
			explicit permutation_matrix(identity_matrix<T>&& mtx) noexcept : identity_matrix<T>(mtx) {}

		public:
			explicit permutation_matrix(
					const square_matrix<T>& mtx,
					const std::size_t row_ind,
					const std::size_t column_ind
			) noexcept : identity_matrix<T>(mtx.rows_number()) {
				this->get_unchecked(row_ind).get_unchecked(row_ind) = 0;
				this->get_unchecked(column_ind).get_unchecked(column_ind) = 0;

				this->get_unchecked(row_ind).get_unchecked(column_ind) = 1;
				this->get_unchecked(column_ind).get_unchecked(row_ind) = 1;
			}

			~permutation_matrix() noexcept = default;

			[[nodiscard]] inline permutation_matrix clone() const noexcept {
				return permutation_matrix(*this);
			}
		};

		template <typename T> inline std::istream& operator >> (std::istream& in, matrix<T>& mtx) noexcept {
			for (auto& cell : mtx)
				in >> cell;

			return in;
		}

		template <typename T> inline std::ostream& operator << (std::ostream& out, const matrix<T>& mtx) noexcept {
			for (auto row_it = mtx.rows_begin(); row_it != mtx.rows_end(); ++row_it) {
				std::copy(row_it->begin(), std::prev(row_it->end()), std::ostream_iterator<T>(out, " "));
				out << *std::prev(row_it->end()) << std::endl;
			}

			return out;
		}

		template <> inline std::ostream& operator << (std::ostream& out, const matrix<double>& mtx) noexcept {
			for (auto row_it = mtx.rows_begin(); row_it != mtx.rows_end(); ++row_it) {
				std::transform(
						row_it->begin(),
						std::prev(row_it->end()),
						std::ostream_iterator<double>(out, " "),
						[](const auto& x) { return x == 0 ? 0 : x; }
				);

				const auto last = *std::prev(row_it->end());
				out << (last == 0 ? 0 : last) << std::endl;
			}

			return out;
		}

		template <typename T> [[nodiscard]] inline bool is_siedel_applicable(const matrix<T>& mtx) noexcept {
			return std::all_of(
					mtx.rows_begin(),
					mtx.rows_end(),
					[i = std::size_t(0)](const auto& row) mutable {
						const auto diagonal_elem = std::abs(row.get_unchecked(i));

						const auto sum = std::accumulate(row.begin(), row.end(), T(0), [&i, q = std::size_t(0)](const auto& acc, const auto& x) mutable {
							return acc + (q++ != i ? std::abs(x) : 0);
						});

						++i;
						return diagonal_elem >= sum;
					}
			);
		}
	}

	namespace column_vec {
		template <typename T> class column_vector : public mtx::matrix<T> {

		public:
			explicit column_vector(const std::size_t size) noexcept : mtx::matrix<T>(size, 1) {}
			explicit column_vector(const mtx::matrix<T>& mtx) noexcept :mtx::matrix<T>(mtx) {}
			explicit column_vector(mtx::matrix<T>&& mtx) noexcept :mtx::matrix<T>(mtx) {}

			column_vector(const std::size_t size, const T& elem) noexcept : mtx::matrix<T>(size, std::vector<T> { elem }) {}
			column_vector(const std::size_t size, T&& elem) noexcept : mtx::matrix<T>(size, std::vector<T> { elem }) {}
			~column_vector() noexcept = default;

			[[nodiscard]] inline std::optional<column_vector> operator+(const column_vector& other) const noexcept {
				const auto rows = this->rows_number();

				if (rows != other.rows_number())
					return std::nullopt;

				column_vector result(rows);

				auto this_it = this->begin();
				auto other_it = other.begin();

				for (auto& elem : result) {
					elem = *this_it + *other_it;
					++this_it, ++other_it;
				}

				return std::make_optional(result);
			}

			[[nodiscard]] inline std::optional<column_vector> operator-(const column_vector& other) const noexcept {
				const auto rows = this->rows_number();

				if (rows != other.rows_number())
					return std::nullopt;

				column_vector result(rows);

				auto this_it = this->begin();
				auto other_it = other.begin();

				for (auto& elem : result) {
					elem = *this_it - *other_it;
					++this_it, ++other_it;
				}

				return std::make_optional(result);
			}

			[[nodiscard]] inline T& get_unchecked(const std::size_t index) noexcept {
				return this->mtx[index].get_unchecked(0);
			}

			[[nodiscard]] inline const T& get_unchecked(const std::size_t index) const noexcept {
				return this->mtx[index].get_unchecked(0);
			}

			[[nodiscard]] inline double norm() const noexcept {
				return std::sqrt(std::accumulate(this->begin(), this->end(), 0.0, [](const auto& acc, const auto& x) {
					return acc + x * x;
				}));
			}
		};
	}
}

template <typename T> inline agla::mtx::matrix<T> read_matrix() noexcept {
	int n = 0, m = 0;
	std::scanf("%d%d", &n, &m);
	agla::mtx::matrix<T> mtx(n, m);
	std::cin >> mtx;
	return mtx;
}

template <typename T> inline agla::mtx::square_matrix<T> read_square_matrix() noexcept {
	int n = 0;
	std::scanf("%d", &n);
	agla::mtx::square_matrix<T> mtx(n);
	std::cin >> mtx;
	return mtx;
}

int main() {
	std::cout << std::setprecision(4) << std::fixed;

	auto a = read_square_matrix<double>();

	std::size_t b_size = 0;
	std::scanf("%zu", &b_size);

	if (!agla::mtx::is_siedel_applicable(a)) {
		std::puts("The method is not applicable!");
		return 0;
	}

	agla::column_vec::column_vector<double> b(b_size);
	std::cin >> b;

	double eps = 0;
	std::scanf("%lf", &eps);

	auto diagonal_a = agla::mtx::square_matrix<double>(a.size(), 0);

	for (std::size_t i = 0; i < a.size(); ++i)
		diagonal_a.get_unchecked(i).get_unchecked(i) = a.get_unchecked(i).get_unchecked(i);

	const auto da_inv_opt = diagonal_a.inversed();

	if (!da_inv_opt.has_value()) {
		std::puts("The method is not applicable!");
		return 0;
	}

	const auto& da_inv = da_inv_opt.value();
	const auto da_inv_mul_a = agla::mtx::square_matrix((da_inv * a).value());
	const auto identity = agla::mtx::identity_matrix<double>(a.size());

	const auto alpha_opt = identity - da_inv_mul_a;

	if (!alpha_opt.has_value()) {
		std::puts("The method is not applicable!");
		return 0;
	}

	auto beta_opt = da_inv * b;

	if (!beta_opt.has_value()) {
		std::puts("The method is not applicable!");
		return 0;
	}

	auto lower_a = a;

	for (std::size_t i = 0; i < a.size(); ++i)
		for (std::size_t q = i + 1; q < a.size(); ++q)
			lower_a.get_unchecked(i).get_unchecked(q) = 0;

	auto lower_a_inv_opt = lower_a.inversed();

	if (!lower_a_inv_opt.has_value()) {
		std::puts("The method is not applicable!");
		return 0;
	}

	auto& lower_a_inv = lower_a_inv_opt.value();

	auto upper_a = a;

	for (std::size_t i = 0; i < a.size(); ++i)
		for (std::size_t q = 0; q <= i; ++q)
			upper_a.get_unchecked(i).get_unchecked(q) = 0;

	auto& alpha = alpha_opt.value();
	auto& beta = beta_opt.value();

	auto lower_a_print = alpha;

	for (std::size_t i = 0; i < a.size(); ++i)
		for (std::size_t q = i + 1; q < a.size(); ++q)
			lower_a_print.get_unchecked(i).get_unchecked(q) = 0;

	auto lower_a_print_inv_opt = lower_a.inversed();

	if (!lower_a_print_inv_opt.has_value()) {
		std::puts("The method is not applicable!");
		return 0;
	}

	auto& lower_a_print_inv = lower_a_inv_opt.value();

	auto upper_a_print = alpha;

	for (std::size_t i = 0; i < a.size(); ++i)
		for (std::size_t q = 0; q <= i; ++q)
			upper_a_print.get_unchecked(i).get_unchecked(q) = 0;

	std::puts("beta:");
	std::cout << beta;

	std::puts("alpha:");
	std::cout << alpha;

	std::puts("B:");
	std::cout << lower_a_print;

	std::puts("C:");
	std::cout << upper_a_print;

	auto i_b = (identity - lower_a_print).value();
	std::puts("I-B:");
	std::cout << i_b;

	auto i_b_inv = i_b.inversed().value();
	std::puts("(I-B)_-1:");
	std::cout << i_b_inv;

	const auto have_to_stop = [&eps](const agla::mtx::matrix<double>& mtx) {
		const auto sum = std::accumulate(mtx.begin(), mtx.end(), 0.0, [](const auto& acc, const auto& x) {
			return acc + x;
		});

		return std::abs(sum - static_cast<double>(mtx.rows_number())) < eps;
	};

	auto x = agla::column_vec::column_vector(std::move(beta));
	std::size_t sol = 0;

	for (;;) {
		std::printf("x(%zu):\n", sol++);
		std::cout << x;

		auto next = agla::column_vec::column_vector(
				(lower_a_inv * (b - agla::column_vec::column_vector(
						(upper_a * x).value()
				)).value()).value()
		);

		const auto cur_eps = (next - x).value().norm();
		std::printf("e: %.4lf\n", cur_eps);
		x = next;

		if (cur_eps <= eps)
			break;
	}

	std::printf("x(%zu):\n", sol);
	std::cout << x;

	return 0;
}
