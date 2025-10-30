"""
Unit tests for datatable utilities.

Test coverage:
- RemoteDataTable filtering
- RemoteDataTable sorting
- RemoteDataTable pagination
"""
import pytest
import pandas as pd

from app.common.utils.datatable_utils import (
    DataTableEvent,
    RemoteDataTable,
    Sort
)


@pytest.mark.unit
@pytest.mark.api
class TestRemoteDataTable:
    """Test RemoteDataTable class."""

    @pytest.fixture
    def sample_dataframe(self):
        """Create sample DataFrame for testing with deep copy to prevent test pollution."""
        df = pd.DataFrame({
            'name': ['Alice', 'Bob', 'Charlie', 'David', 'Eve'],
            'age': [25, 30, 35, 28, 32],
            'city': ['New York', 'Los Angeles', 'Chicago', 'Houston', 'Phoenix'],
            'score': [85, 90, 75, 88, 92]
        })
        return df.copy(deep=True)

    @pytest.fixture
    def sample_event(self):
        """Create sample DataTableEvent for testing."""
        return DataTableEvent(
            file_name="test.csv",
            page=1,
            perPage=10,
            columnFilters={},
            sort=[]
        )

    def test_remote_datatable_initialization(self, sample_dataframe, sample_event):
        """Test RemoteDataTable initialization."""
        # Act
        remote_table = RemoteDataTable(sample_dataframe, sample_event)

        # Assert
        assert len(remote_table.df) == 5
        assert remote_table.total_records == 5
        assert remote_table.server_params == sample_event
        assert remote_table.filtered_sorted_df.equals(sample_dataframe)

    def test_apply_filters_single_column(self, sample_dataframe, sample_event):
        """Test filtering with single column filter."""
        # Arrange
        sample_event.columnFilters = {'name': 'Alice'}
        remote_table = RemoteDataTable(sample_dataframe, sample_event)

        # Act
        remote_table.apply_filters()

        # Assert
        assert len(remote_table.filtered_sorted_df) == 1
        assert remote_table.filtered_sorted_df.iloc[0]['name'] == 'Alice'

    def test_apply_filters_multiple_columns(self, sample_dataframe, sample_event):
        """Test filtering with multiple column filters."""
        # Arrange
        sample_event.columnFilters = {'city': 'New', 'name': 'Alice'}
        remote_table = RemoteDataTable(sample_dataframe, sample_event)

        # Act
        remote_table.apply_filters()

        # Assert
        assert len(remote_table.filtered_sorted_df) == 1
        assert remote_table.filtered_sorted_df.iloc[0]['name'] == 'Alice'
        assert remote_table.filtered_sorted_df.iloc[0]['city'] == 'New York'

    def test_apply_filters_no_match(self, sample_dataframe, sample_event):
        """Test filtering with no matching results."""
        # Arrange
        sample_event.columnFilters = {'name': 'NonexistentName'}
        remote_table = RemoteDataTable(sample_dataframe, sample_event)

        # Act
        remote_table.apply_filters()

        # Assert
        assert len(remote_table.filtered_sorted_df) == 0

    def test_apply_sorting_ascending(self, sample_dataframe, sample_event):
        """Test sorting in ascending order."""
        # Arrange
        sample_event.sort = [Sort(field='age', type='asc')]
        remote_table = RemoteDataTable(sample_dataframe, sample_event)

        # Act
        remote_table.apply_sorting()

        # Assert
        ages = remote_table.filtered_sorted_df['age'].tolist()
        assert ages == [25, 28, 30, 32, 35]

    def test_apply_sorting_descending(self, sample_dataframe, sample_event):
        """Test sorting in descending order."""
        # Arrange
        sample_event.sort = [Sort(field='score', type='desc')]
        remote_table = RemoteDataTable(sample_dataframe, sample_event)

        # Act
        remote_table.apply_sorting()

        # Assert
        scores = remote_table.filtered_sorted_df['score'].tolist()
        assert scores == [92, 90, 88, 85, 75]

    def test_apply_sorting_no_sort(self, sample_dataframe, sample_event):
        """Test behavior when no sorting is specified."""
        # Arrange
        sample_event.sort = []
        remote_table = RemoteDataTable(sample_dataframe, sample_event)

        # Act
        remote_table.apply_sorting()

        # Assert - data should remain in original order
        assert remote_table.filtered_sorted_df.equals(sample_dataframe)

    def test_get_paginated_data_first_page(self, sample_dataframe, sample_event):
        """Test pagination for first page."""
        # Arrange
        sample_event.page = 1
        sample_event.perPage = 2
        remote_table = RemoteDataTable(sample_dataframe, sample_event)

        # Act
        result = remote_table.get_paginated_data()

        # Assert
        assert result['totalRecords'] == 5
        assert len(result['df']) == 2
        assert result['df'].iloc[0]['name'] == 'Alice'
        assert result['df'].iloc[1]['name'] == 'Bob'

    def test_get_paginated_data_second_page(self, sample_dataframe, sample_event):
        """Test pagination for second page."""
        # Arrange
        sample_event.page = 2
        sample_event.perPage = 2
        remote_table = RemoteDataTable(sample_dataframe, sample_event)

        # Act
        result = remote_table.get_paginated_data()

        # Assert
        assert result['totalRecords'] == 5
        assert len(result['df']) == 2
        assert result['df'].iloc[0]['name'] == 'Charlie'

    def test_get_paginated_data_last_page_partial(self, sample_dataframe, sample_event):
        """Test pagination for last page with partial results."""
        # Arrange
        sample_event.page = 3
        sample_event.perPage = 2
        remote_table = RemoteDataTable(sample_dataframe, sample_event)

        # Act
        result = remote_table.get_paginated_data()

        # Assert
        assert result['totalRecords'] == 5
        assert len(result['df']) == 1  # Only one record on last page
        assert result['df'].iloc[0]['name'] == 'Eve'

    def test_get_filtered_sorted_paginated_data_combined(self, sample_dataframe, sample_event):
        """Test combined filtering, sorting, and pagination."""
        # Arrange
        sample_event.columnFilters = {'city': 'Angeles'}  # Matches only "Los Angeles"
        sample_event.sort = [Sort(field='age', type='asc')]
        sample_event.page = 1
        sample_event.perPage = 2
        remote_table = RemoteDataTable(sample_dataframe, sample_event)

        # Act
        result = remote_table.get_filtered_sorted_paginated_data()

        # Assert
        assert result['totalRecords'] == 5  # Original total
        assert len(result['df']) == 1  # Only one match after filtering
        # Should be sorted by age ascending after filtering
        ages = result['df']['age'].tolist()
        assert ages[0] == 30  # Bob (Los Angeles)

    def test_apply_filters_numeric_column(self, sample_dataframe, sample_event):
        """Test filtering on numeric column (converted to string for contains check)."""
        # Arrange
        sample_event.columnFilters = {'age': '3'}  # Matches ages containing '3' (30, 35, 32)
        remote_table = RemoteDataTable(sample_dataframe, sample_event)

        # Act
        remote_table.apply_filters()

        # Assert
        assert len(remote_table.filtered_sorted_df) == 3
        ages = set(remote_table.filtered_sorted_df['age'].tolist())
        assert ages == {30, 32, 35}

    def test_apply_filters_case_sensitivity(self, sample_dataframe, sample_event):
        """Test case sensitivity in filtering (pandas str.contains is case-sensitive by default)."""
        # Arrange
        sample_event.columnFilters = {'name': 'alice'}  # lowercase 'alice'
        remote_table = RemoteDataTable(sample_dataframe, sample_event)

        # Act
        remote_table.apply_filters()

        # Assert - Should not match 'Alice' (case-sensitive)
        assert len(remote_table.filtered_sorted_df) == 0

    def test_apply_filters_case_insensitive_partial_match(self, sample_dataframe, sample_event):
        """Test partial match filtering works correctly."""
        # Arrange
        sample_event.columnFilters = {'city': 'hoe'}  # Matches 'Phoenix' (contains 'hoe')
        remote_table = RemoteDataTable(sample_dataframe, sample_event)

        # Act
        remote_table.apply_filters()

        # Assert
        assert len(remote_table.filtered_sorted_df) == 1
        assert remote_table.filtered_sorted_df.iloc[0]['city'] == 'Phoenix'

    def test_apply_sorting_multi_column_only_uses_first(self, sample_dataframe, sample_event):
        """Test that multi-column sort only uses the first sort specification."""
        # Arrange - Add multiple sorts, but implementation only uses first
        sample_event.sort = [
            Sort(field='age', type='desc'),  # This should be applied
            Sort(field='name', type='asc')   # This should be ignored
        ]
        remote_table = RemoteDataTable(sample_dataframe, sample_event)

        # Act
        remote_table.apply_sorting()

        # Assert - Only sorted by age descending
        ages = remote_table.filtered_sorted_df['age'].tolist()
        assert ages == [35, 32, 30, 28, 25]
        # Names are not sorted alphabetically within same age groups

    def test_get_paginated_data_page_beyond_data(self, sample_dataframe, sample_event):
        """Test pagination when requesting page beyond available data."""
        # Arrange
        sample_event.page = 10  # Way beyond available data
        sample_event.perPage = 2
        remote_table = RemoteDataTable(sample_dataframe, sample_event)

        # Act
        result = remote_table.get_paginated_data()

        # Assert - Should return empty DataFrame
        assert result['totalRecords'] == 5
        assert len(result['df']) == 0

    def test_get_paginated_data_large_perpage(self, sample_dataframe, sample_event):
        """Test pagination with perPage larger than total records."""
        # Arrange
        sample_event.page = 1
        sample_event.perPage = 100  # Larger than dataset
        remote_table = RemoteDataTable(sample_dataframe, sample_event)

        # Act
        result = remote_table.get_paginated_data()

        # Assert - Should return all records
        assert result['totalRecords'] == 5
        assert len(result['df']) == 5

    def test_total_records_unchanged_after_filtering(self, sample_dataframe, sample_event):
        """Test that totalRecords reflects original dataset size, not filtered size."""
        # Arrange
        sample_event.columnFilters = {'name': 'Alice'}  # Filters to 1 record
        remote_table = RemoteDataTable(sample_dataframe, sample_event)

        # Act
        result = remote_table.get_filtered_sorted_paginated_data()

        # Assert - totalRecords should still be 5 (original), not 1 (filtered)
        assert result['totalRecords'] == 5
        assert len(result['df']) == 1  # But only 1 record in paginated data

    # ===== EDGE CASE TESTS =====

    def test_get_paginated_data_zero_perpage(self, sample_dataframe, sample_event):
        """Test pagination with perPage=0 returns empty DataFrame.

        Edge case: When perPage=0, pagination should return empty result.
        Calculation: start=(1-1)*0=0, end=0+0=0 → iloc[0:0] → empty DataFrame

        This documents pandas iloc behavior where start==end returns empty slice.
        """
        # Arrange
        sample_event.page = 1
        sample_event.perPage = 0
        remote_table = RemoteDataTable(sample_dataframe, sample_event)

        # Act
        result = remote_table.get_paginated_data()

        # Assert
        assert result['totalRecords'] == 5  # Original size unchanged
        assert len(result['df']) == 0  # Empty page due to perPage=0
        assert isinstance(result['df'], pd.DataFrame)  # Still returns DataFrame

    def test_get_paginated_data_negative_perpage(self, sample_dataframe, sample_event):
        """Test pagination with negative perPage value behavior.

        Edge case: perPage=-1 creates negative slicing in pandas iloc.
        Calculation: start=(1-1)*(-1)=0, end=0+(-1)=-1 → iloc[0:-1] → all but last row

        This documents pandas iloc behavior with negative indices.
        Result: Returns rows 0 through -1 (all rows except the last one).
        """
        # Arrange
        sample_event.page = 1
        sample_event.perPage = -1
        remote_table = RemoteDataTable(sample_dataframe, sample_event)

        # Act
        result = remote_table.get_paginated_data()

        # Assert
        assert result['totalRecords'] == 5  # Original size unchanged
        # iloc[0:-1] returns all rows except the last one
        assert len(result['df']) == 4  # All rows except Eve
        assert 'Eve' not in result['df']['name'].tolist()

    def test_get_paginated_data_page_zero(self, sample_dataframe, sample_event):
        """Test pagination with page=0 (invalid page number).

        Edge case: page=0 creates negative start index.
        Calculation: start=(0-1)*10=-10, end=-10+10=0 → iloc[-10:0] → empty

        This documents pandas iloc behavior with negative start index.
        Result: iloc[-10:0] with 5-row DataFrame returns empty (no rows before index 0).
        """
        # Arrange
        sample_event.page = 0
        sample_event.perPage = 10
        remote_table = RemoteDataTable(sample_dataframe, sample_event)

        # Act
        result = remote_table.get_paginated_data()

        # Assert
        assert result['totalRecords'] == 5  # Original size unchanged
        # iloc[-10:0] returns empty for 5-row DataFrame
        assert len(result['df']) == 0
        assert isinstance(result['df'], pd.DataFrame)

    def test_apply_sorting_invalid_direction(self, sample_dataframe, sample_event):
        """Test sorting with invalid sort type (not 'asc' or 'desc').

        Edge case: sort_type='invalid' → ascending=False (default falsy behavior)
        Implementation: ascending = (sort_type == 'asc')
        Result: Any non-'asc' value including 'invalid' is treated as descending.

        This documents that implementation treats unknown sort types as descending.
        """
        # Arrange
        sample_event.sort = [Sort(field='age', type='invalid')]
        remote_table = RemoteDataTable(sample_dataframe, sample_event)

        # Act
        remote_table.apply_sorting()

        # Assert - Should not crash, treats 'invalid' as descending
        ages = remote_table.filtered_sorted_df['age'].tolist()
        assert ages == [35, 32, 30, 28, 25]  # Descending order (ascending=False)

    def test_apply_sorting_nonexistent_field(self, sample_dataframe, sample_event):
        """Test sorting with field that doesn't exist in DataFrame.

        Edge case: Pandas sort_values() raises KeyError for missing columns.
        This test documents expected error behavior for invalid field names.

        Result: KeyError is raised, no silent failure or data corruption.
        """
        # Arrange
        sample_event.sort = [Sort(field='nonexistent_column', type='asc')]
        remote_table = RemoteDataTable(sample_dataframe, sample_event)

        # Act & Assert - Should raise KeyError
        with pytest.raises(KeyError):
            remote_table.apply_sorting()

    def test_apply_filters_empty_value(self, sample_dataframe, sample_event):
        """Test filtering with empty string filter value.

        Edge case: Empty string '' in str.contains() matches all rows.
        Pandas behavior: df[column].str.contains('') returns all True.

        This documents pandas str.contains('') behavior - empty pattern matches everything.
        Result: All rows pass the filter when filter value is empty string.
        """
        # Arrange
        sample_event.columnFilters = {'name': ''}  # Empty filter value
        remote_table = RemoteDataTable(sample_dataframe, sample_event)

        # Act
        remote_table.apply_filters()

        # Assert - Empty string matches all rows
        assert len(remote_table.filtered_sorted_df) == 5  # All rows returned
        assert set(remote_table.filtered_sorted_df['name']) == {'Alice', 'Bob', 'Charlie', 'David', 'Eve'}
