#ifndef BASIC2DMATRIX_H
#define BASIC2DMATRIX_H

template <class Object>
class Basic2DMatrix
{
public:
	Basic2DMatrix( int rows, int cols ) : array( rows )
	{
		for( int i = 0; i < rows; i++ )
			array[ i ].resize( cols );
	}

	Basic2DMatrix( const Basic2DMatrix & rhs ) : array( rhs.array ) { }

	const std::vector<Object> & operator[]( int row ) const
	{ return array[ row ]; }

	std::vector<Object> & operator[]( int row )
	{ return array[ row ]; }

	int numrows( ) const
	{ return array.size( ); }

	int numcols( ) const
	{ return numrows( ) ? array[ 0 ].size( ) : 0; }

private:
	std::vector< std::vector<Object> > array;
};

#endif  // BASIC2DMATRIX

