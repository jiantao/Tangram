#!/usr/bin/python
import sys
import struct
import argparse
import simplejson as json
from pprint import pprint
import cStringIO


"""
Since the lib_table.dat and hist.dat generated by tangram_scan is in binary format,
it is quite hard to check or change the information in them. This program provides
the function to read, print or write the information in these two files. To use it as
a stand alone program, run 'python tangram_view_scan_file.py -h' to see the usage information.
To use it as API to change or rewrite the information in the two files, import the reader and writer
classes.

What should we pay attention to in the lib_table?
1) Usually the special reference name (special: {'name':[...]}). If you do not see them, MEI detection will fail.
2) Crappy read groups. The frag_len_info dict store the stats for each fragment length distribution from each read group.
   If, for some reason, you see the frag_max value is very huge, say > 10000, this read group might be low-quality.

Why do we want to change the information in lib_table.dat?
1) Manually filtered out a certain read group: set the frag_max, frag_med and frag_min to zero.
2) Change the cutoff and trim rate of a read group (this requires loading hist.dat and recalculate the percentile).

To be used by trangram_merge or detect, the output format of lib_table or hist object must be binary.
"""


class LibTableReader(object):
	'''
	Provid the function to read the information in the lib_table.dat file.
	The information will stored in a dict like this:
	{'anchor': {'length': [249250621,
						   -1],
				'md5': ['1b22b98cdeb4a9304cb5d48026a85128',
						'5b6a4b3a81a2d3c134b7d14bf6ad39f1'],
				'name': ['1',
						 'hs37d5'],
				'special': {'name': ['L1', 'AL', 'HE', 'SV', 'PO']},
				'special_prefix': ''},
	 'read_grp': {'cutoff': 0.01,
				  'frag_len_info': [{'frag_max': 664,
									 'frag_med': 499,
									 'frag_min': 326}],
				  'frag_len_max': 664,
				  'name': ['ZIIAREPXRA8'],
				  'sample_map': [0],
				  'seq_tech': [0],
				  'trim_rate': 0.002},
	 'sample': {'name': ['unknown']}}

	'anchor' stores the information of references.
	'sample' stores the information of samples.
	'read_grp' stores the information of read groups.

	'''
	def __init__(self):
		self.lib_table_input = None
		self.anchor_size = 0
		self.sample_size = 0
		self.read_grp_size = 0

	def parse_lib_table_file(self, lib_table_file):
		self.lib_table_input = open(lib_table_file, 'rb')
		# read the number of references, samples and read groups
		self.anchor_size = struct.unpack('<I', self.lib_table_input.read(4))[0]
		self.sample_size = struct.unpack('<I', self.lib_table_input.read(4))[0]
		self.read_grp_size = struct.unpack('<I', self.lib_table_input.read(4))[0]

		lib_table = {}

		lib_table['anchor'] = self._read_anchors()
		lib_table['sample'] = self._read_samples()
		lib_table['read_grp'] = self._read_read_grps()
		lib_table['anchor']['special'] = self._read_specials()

		self.lib_table_input.close()

		return lib_table

	def _read_anchors(self):
		"""
		Reading regular chromosome (not including MEI reference)
		information from the lib_table.dat file.
		"""
		anchor_info = {}
		# length of each chr
		# if the length is negative then those chr won't
		# be used for detection
		block = self.lib_table_input.read(4 * self.anchor_size)
		anchor_info['length'] = list(struct.unpack('<' + 'i' * self.anchor_size, block))

		# md5 string of each chr
		block = self.lib_table_input.read(32 * self.anchor_size)
		anchor_info['md5'] = []
		for i in xrange(0, 32 * self.anchor_size, 32):
			anchor_info['md5'].append(struct.unpack('32s', block[i:i+32])[0])

		# name of each chr
		anchor_info['name'] = []
		for i in xrange(self.anchor_size):
			name_len = struct.unpack('<I', self.lib_table_input.read(4))[0]
			anchor_name = struct.unpack(str(name_len) + 's', self.lib_table_input.read(name_len))[0]
			anchor_info['name'].append(anchor_name)

		# prefix of special chr (MEI) name
		prefix_len = struct.unpack('<I', self.lib_table_input.read(4))[0]
		anchor_info['special_prefix'] = ''
		if prefix_len > 0:
			anchor_info['special_prefix'] = struct.unpack(str(prefix_len) + 's', self.lib_table_input.read(prefix_len))[0]

		return anchor_info

	def _read_samples(self):
		"""
		Reading sample information from the lib_table.dat file.
		"""
		sample_info = {}

		# name of each sample
		sample_info['name'] = []
		for i in xrange(self.sample_size):
			name_len = struct.unpack('I', self.lib_table_input.read(4))[0]
			sample_name = struct.unpack(str(name_len) + 's', self.lib_table_input.read(name_len))[0]
			sample_info['name'].append(sample_name)

		return sample_info

	def _read_read_grps(self):
		"""
		Reading read group information from the lib_table.dat file.
		"""
		read_grp_info = {}

		# map read group to its corresponding sample
		block = self.lib_table_input.read(4 * self.read_grp_size)
		read_grp_info['sample_map'] = list(struct.unpack('<' + str(self.read_grp_size) + 'i', block))

		# name of each read group
		read_grp_info['name'] = []
		for i in xrange(self.read_grp_size):
			name_len = struct.unpack('I', self.lib_table_input.read(4))[0]
			read_grp_name = struct.unpack(str(name_len) + 's', self.lib_table_input.read(name_len))[0]
			read_grp_info['name'].append(read_grp_name)

		# largest fragment length among all the read groups
		# if this value is zero then there are nothing to analyze
		# and we will bug out here
		read_grp_info['frag_len_max'] = struct.unpack('<I', self.lib_table_input.read(4))[0]
		if read_grp_info['frag_len_max'] == 0:
			return read_grp_info

		# the cutoff percentage of the histogram of fragment length of normal
		# read pairs.
		read_grp_info['cutoff'] = struct.unpack('<d', self.lib_table_input.read(8))[0]
		# how much percent of data are pre-trimmed in fragnment length
		# distribution to exclude outliers
		read_grp_info['trim_rate'] = struct.unpack('<d', self.lib_table_input.read(8))[0]

		# sequencing technology of each read group
		block = self.lib_table_input.read(self.read_grp_size)
		read_grp_info['seq_tech'] = list(struct.unpack(str(self.read_grp_size) + 'b', block))

		# for each read group get its upper limit, lower limit (defined by
		# cutoff) and median of the fragment length distribution
		read_grp_info['frag_len_info'] = []
		for i in xrange(self.read_grp_size):
			frag_len_info = struct.unpack('<3i', self.lib_table_input.read(4*3))
			frag_len_dict = {}
			frag_len_dict['frag_med'] = frag_len_info[0]
			frag_len_dict['frag_max'] = frag_len_info[1]
			frag_len_dict['frag_min'] = frag_len_info[2]
			read_grp_info['frag_len_info'].append(frag_len_dict)

		return read_grp_info

	def _read_specials(self):
		"""
		Reading special reference information from the lib_table.dat file.
		"""
		special_info = {}

		# number of special references
		special_size = struct.unpack('<I', self.lib_table_input.read(4))[0]
		if special_size == 0:
			return special_info

		# name of each special reference
		special_info['name'] = []
		for i in xrange(special_size):
			special_name = struct.unpack('2s', self.lib_table_input.read(2))[0]
			special_info['name'].append(special_name)

		return special_info

class LibTableWriter(object):
	def __init__(self):
		self.bin_output = None

	def write_lib_table(self, lib_table, out_type, output):
		"""
		Write the lib table to the output flow in text, json or binary format.
		"""
		if out_type == 'text':
			pprint(lib_table, output)
		elif out_type == 'json':
			output.write(json.dumps(lib_table))
		elif out_type == 'binary':
			self.bin_output = cStringIO.StringIO()
			self._generate_binary_output(lib_table)
			output.write(self.bin_output.getvalue())
			self.bin_output.close()

		if output != sys.stdout:
			output.close()

	def _generate_binary_output(self, lib_table):
		"""
		Generate the binary output string from a lib table object.
		"""
		anchor_size = len(lib_table['anchor']['name'])
		sample_size = len(lib_table['sample']['name'])
		read_grp_size = len(lib_table['read_grp']['name'])
		if anchor_size == 0 or sample_size == 0 or read_grp_size == 0:
			raise ValueError('Unable to write an empty lib table.')

		self.bin_output.write(struct.pack('<I', anchor_size))
		self.bin_output.write(struct.pack('<I', sample_size))
		self.bin_output.write(struct.pack('<I', read_grp_size))

		self._generate_anchor_binary_output(lib_table['anchor'])
		self._generate_sample_binary_output(lib_table['sample'])
		self._generate_read_grp_binary_output(lib_table['read_grp'])
		self._generate_special_binary_output(lib_table['anchor']['special'])

	def _generate_anchor_binary_output(self, anchor_info):
		"""
		Convert the anchor information into binary output.
		"""
		anchor_size = len(anchor_info['name'])
		# check the data consistency of name, length and md5 of references
		if anchor_size != len(anchor_info['md5']) or anchor_size != len(anchor_info['length']):
			raise ValueError('The number of items for anchor length, md5 strings and names are not consistent with each other.')

		# write the reference length
		self.bin_output.write(struct.pack('<' + str(anchor_size) + 'i', *anchor_info['length']))
		# write the reference md5 string
		self.bin_output.write(struct.pack('32s' * anchor_size, *anchor_info['md5']))
		# write the reference name
		self._convert_name_list(anchor_info['name'])

		# write the prefix string of special reference
		special_prefix_len = len(anchor_info['special_prefix'])
		self.bin_output.write(struct.pack('<I', special_prefix_len))
		if special_prefix_len > 0:
			self.bin_output.write(struct.pack(str(special_prefix_len) + 'I', anchor_info['special_prefix']))

	def _generate_sample_binary_output(self, sample_info):
		"""
		Convert the sample information into binary output.
		"""
		self._convert_name_list(sample_info['name'])

	def _generate_read_grp_binary_output(self, read_grp_info):
		"""
		Convert the read group information into binary output.
		"""
		# the length of all lists, name, frag_len_info and seq_tech, should
		# have the same length
		read_grp_size = len(read_grp_info['name'])
		if read_grp_size != len(read_grp_info['sample_map']) or \
		   read_grp_size != len(read_grp_info['frag_len_info']) or \
		   read_grp_size != len(read_grp_info['seq_tech']):
			raise ValueError('The number of items for read group sample map, fragment length information and names are not consistent with each other.')

		# write the sample map (map the read group to its corrponding sample)
		self.bin_output.write(struct.pack('<' + str(read_grp_size) + 'i', *read_grp_info['sample_map']))
		# write the name of each read group
		self._convert_name_list(read_grp_info['name'])
		# write the largest fragment length among all read groups
		self.bin_output.write(struct.pack('<I', read_grp_info['frag_len_max']))
		# the frag_len_max value is zero we will skip the rest of information
		# in the read_grp_info project
		if read_grp_info['frag_len_max'] == 0:
			return

		# write the fragment length distribution cutoff
		self.bin_output.write(struct.pack('<d', read_grp_info['cutoff']))
		# write the fragment length distribution trim rate
		self.bin_output.write(struct.pack('<d', read_grp_info['trim_rate']))
		# write the sequencing technology for each read group
		self.bin_output.write(struct.pack(str(read_grp_size) + 'b', *read_grp_info['seq_tech']))

		# write the fragment length distribution stats
		for frag_len_dict in read_grp_info['frag_len_info']:
			self.bin_output.write(struct.pack('<i', frag_len_dict['frag_med']))
			self.bin_output.write(struct.pack('<i', frag_len_dict['frag_max']))
			self.bin_output.write(struct.pack('<i', frag_len_dict['frag_min']))

	def _generate_special_binary_output(self, special_info):
		"""
		Convert the special reference info into binary output.
		"""
		special_size = len(special_info['name'])
		self.bin_output.write(struct.pack('<I', special_size))
		if special_size == 0:
			return

		self.bin_output.write(struct.pack('2s' * special_size, *special_info['name']))

	def _convert_name_list(self, name_list):
		"""
		Convert a list of names (strings) into binary output.
		"""
		for name in name_list:
			name_len = len(name)
			# name should not be empty
			if name_len == 0:
				raise ValueError("Invalid sample name.")
			self.bin_output.write(struct.pack('<I', name_len))
			self.bin_output.write(struct.pack(str(name_len) + 's', name))


class HistReader(object):
	"""
	Provide the function of reading the information in the hist.dat file.
	The information in hist.dat looks like:
	[{'bin': [107,
			  50091896],
	  'freq': [1,
			   1]}]

	This is a list of histograms (dict{'bin':[...], 'freq':[...]}). Each histogram
	is generated from each read group. The order of histograms is the same as the
	read groups in the lib_table object.
	"""
	def __init__(self):
		self.hist_input = None
		self.hist_size = 0

	def parse_hist(self, hist_file):
		"""
		Read the hist.dat file generated by tangram_scan.
		"""
		self.hist_input = open(hist_file, 'rb')
		self.hist_size = struct.unpack('<I', self.hist_input.read(4))[0]
		hist_table = self._read_hist_array()

		return hist_table

	def _read_hist_array(self):
		"""
		Read each hist array from the binary file.
		"""
		hist_array = []
		for i in xrange(self.hist_size):
			hist_len = struct.unpack('<I', self.hist_input.read(4))[0]
			if hist_len == 0:
				hist_array.append({'bin':None, 'freq':None})
				continue

			hist_element = {}
			hist_element['bin'] = list(struct.unpack(str(hist_len) + 'I', self.hist_input.read(4*hist_len)))
			hist_element['freq'] = list(struct.unpack(str(hist_len) + 'Q', self.hist_input.read(8*hist_len)))
			hist_array.append(hist_element)

		return hist_array


class HistTableWriter(object):
	def __init__(self):
		self.bin_output = None

	def write_hist_table(self, hist_table, out_type, output):
		"""
		Write the hist table to the output flow in text, json or binary format.
		"""
		if out_type == 'text':
			pprint(hist_table, output)
		elif out_type == 'json':
			output.write(json.dumps(hist_table))
		elif out_type == 'binary':
			self.bin_output = cStringIO.StringIO()
			self._generate_binary_output(hist_table)
			output.write(self.bin_output.getvalue())
			self.bin_output.close()

		if output != sys.stdout:
			output.close()

	def _generate_binary_output(self, hist_table):
		"""
		Generate the binary output string from a hist table object.
		"""
		hist_size = len(hist_table)
		if hist_size == 0:
			raise ValueError('Unable to write an empty hist table.')

		self.bin_output.write(struct.pack('<I', hist_size))
		for i in xrange(hist_size):
			if hist_table[i]['bin'] is None:
				self.bin_output.write(struct.pack('<I', 0))
				continue

			hist_len = len(hist_table[i]['bin'])
			if hist_len != len(hist_table[i]['freq']):
				raise ValueError('The number of bin is not consistent with the number of frequency in a histogram.')

			self.bin_output.write(struct.pack('<I', hist_len))
			self.bin_output.write(struct.pack('<' + str(hist_len) + 'I', *hist_table[i]['bin']))
			self.bin_output.write(struct.pack('<' + str(hist_len) + 'Q', *hist_table[i]['freq']))


def get_cmd_args():
	parser = argparse.ArgumentParser(description='Showing the content in the binary files generated by tangram_scan.')
	group = parser.add_mutually_exclusive_group(required=True)
	group.add_argument('--lib', '-l', help='input lib_table.dat file')
	group.add_argument('--hist', '-s', help='input hist.dat file')
	parser.add_argument('--out', '-o', default=sys.stdout, type=argparse.FileType('w'), help='output file name [default stdin]')
	parser.add_argument('--type', '-t', default='text', choices=['text', 'json', 'binary'],
			help='choose on of the following output fomat: text, json or binary [default: text]')

	args = parser.parse_args()
	return args


if __name__ == '__main__':
	args = get_cmd_args()

	if args.lib is not None:
		lib_table_reader = LibTableReader()
		lib_table = lib_table_reader.parse_lib_table_file(args.lib)
		lib_table_writer = LibTableWriter()
		lib_table_writer.write_lib_table(lib_table, args.type, args.out)

	if args.hist is not None:
		hist_reader = HistReader()
		hist_table = hist_reader.parse_hist(args.hist)
		hist_table_writer = HistTableWriter()
		hist_table_writer.write_hist_table(hist_table, args.type, args.out)
