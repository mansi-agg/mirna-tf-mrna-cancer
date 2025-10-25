import tabula
import pandas as pd

# # Extract tables from a PDF
# tables = tabula.read_pdf("data_files/Cellosaurus query breast cancer.pdf", pages='all')
# # tables = tabula.read_pdf("data_files/Cellosaurus query melanoma.pdf", pages='all')

# # # Print the first table
# # print(tables[81])

# cleaned_tables = []

# for table in tables:
#     # Create an empty DataFrame to store the cleaned table
#     cleaned_table = []

#     header_as_row = list(table.columns)
#     cleaned_table.append(header_as_row)

#     i = 0
#     while i < len(table):
#         # Check if the current row is a "split row" (next row has values that should be merged)
#         if pd.isna(table.iloc[i, 1]) and i + 1 < len(table):
#             # Merge the current row with the next row
#             new_row = [
#                 table.iloc[i, 0],  # First column stays the same
#                 table.iloc[i + 1, 1],  # Take the second column from the next row
#                 table.iloc[i + 1, 2]  # Take the third column from the next row
#             ]
#             cleaned_table.append(new_row)
#             i += 2  # Skip the next row since it's merged
#         else:
#             # If it's a regular row, just take it as is
#             cleaned_table.append(table.iloc[i].values)
#             i += 1

#     # Convert the cleaned table into a DataFrame and add it to the list
#     cleaned_df = pd.DataFrame(cleaned_table)

#     cleaned_tables.append(cleaned_df)

# # Combine all cleaned tables into a single DataFrame
# final_df = pd.concat(cleaned_tables, ignore_index=True)
# final_df = final_df.reset_index(drop=True)
# final_df = final_df.drop_duplicates()
# final_df.columns = ['Accession', 'Cell line name', 'Species of origin']

# final_df.to_csv("data_files/Cellosaurus query breast cancer.csv", index_label=False)