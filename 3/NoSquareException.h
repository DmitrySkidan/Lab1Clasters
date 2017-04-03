#pragma once

#include <string>
#include <stdexcept>

namespace matrix
{

	class NoSquareException : public std::exception
	{

	public:
		NoSquareException();

		NoSquareException(const std::wstring &message);

	};

}
