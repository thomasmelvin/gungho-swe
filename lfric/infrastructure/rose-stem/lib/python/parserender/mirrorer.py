##############################################################################
# Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
# For further details please refer to the file LICENCE.original which you
# should have received as part of this distribution.
##############################################################################
"""
Upload a filetree to a server of some sort with arbitrary fileTransforms
en-route.
"""
from abc import ABC, abstractmethod
import ftplib
import io
import netrc
import os
import os.path
import re
import shutil
import stat
from urllib.parse import urlparse
import xml.etree.ElementTree as ET


##############################################################################
class TreeVisitor(ABC):  # pylint: disable=too-few-public-methods
    """
    Interface for file tree visitors.
    """
    @abstractmethod
    def visit(self, directory, subdirs, files):
        """
        Examines a directory.
        """
        raise NotImplementedError()


##############################################################################
class GenerateTreeIndices(TreeVisitor):
    """
    Create index files where needed.
    """
    def __init__(self, renderer):
        self._renderer = renderer
        self._directory_list = {}
        self._index_found = -999

    def visit(self, directory, subdirs, files):
        tree_level = directory.count(os.sep)
        if tree_level < self._index_found:
            self._index_found = -999

        if 'index.html' in files:
            self._index_found = tree_level
            return

        if not self._index_found:
            self._directory_list[directory] = subdirs + files

    def new_files(self):
        """
        It's too hot to work out what this does.
        """
        page_list = {}
        for filename, file_list in self._directory_list.items():
            content = self._renderer.render('"{}" Directory'.format(filename),
                                            file_list)
            page_list[os.path.join(filename, 'index.html')] = content
        return page_list


##############################################################################
class XmlTransformation(ABC):  # pylint: disable=too-few-public-methods
    """
    Interface for transforming XML documents.
    """
    @abstractmethod
    def transform(self, root):
        """
        Transform the XML element tree in some way.
        """
        raise NotImplementedError()


##############################################################################
class RemoveTagsWithClass(XmlTransformation):
    """
    Removes any tag with the specified class from an XHTML document.
    """
    def __init__(self, class_name):
        self._class_name = class_name

    def transform(self, root):
        # TODO: This does not handle the case where there is
        #       more than one class specified in the "class" attribute.
        xpath = ".//*[@class='{}']".format(self._class_name)
        for parent in root.findall(xpath + '/..'):
            for element in parent.findall(xpath):
                parent.remove(element)


##############################################################################
class AddRelativePrefix(XmlTransformation):
    """
    Adds a prefix to all relative URLs found in HREF attributes.
    """
    def __init__(self, prefix):
        self._prefix = prefix
        if self._prefix[-1] != '/':
            self._prefix = self._prefix + '/'

    def transform(self, root):
        xpath = './/*[@href]'
        for parent in root.findall(xpath + '/..'):
            for element in parent.findall(xpath):
                url = urlparse(element.attrib['href'])
                if not url.scheme and not url.netloc:
                    with_prefix = self._prefix + element.attrib['href']
                    element.attrib['href'] = with_prefix


##############################################################################
class NoNakedURLs(XmlTransformation):  # pylint: disable=too-few-public-methods
    """
    Any URLs which do not specify a file have a default index page added.
    This assumes that all files have an extension.
    """
    def __init__(self, page_name):
        self._default_page_name = page_name

    def transform(self, root):
        xpath = r'.//*[@href]'
        file_pattern = re.compile(r'.*\..+$')
        for parent in root.findall(xpath + '/..'):
            for element in parent.findall(xpath):
                if not file_pattern.match(element.attrib['href']):
                    if not element.attrib['href'].endswith('/'):
                        element.attrib['href'] = element.attrib['href'] + '/'
                    element.attrib['href'] += 'index.html'


##############################################################################
class Credentials(ABC):  # pylint: disable=too-few-public-methods
    """
    Interface for fetching access credentials.
    """
    @abstractmethod
    def get_credentials(self):
        """
        Gets username and password.
        """
        raise NotImplementedError()


##############################################################################
class NetrcCredentials(Credentials):  # pylint: disable=too-few-public-methods
    """
    Extracts credentials from a .netrc file.
    """
    def __init__(self, host):
        """
        Construct a Credentials object which will use a .netrc file as its
        source.

        host - Name of the host as it appears in .netrc
        """
        self._host = host

    def get_credentials(self):
        """
        Returns the first match found.

        The file is opened and read here rather than in the constructor to
        minimise the amount of time key matter spends in memory. Hopefully.
        This whole thing is a bit shoddy.
        """
        credentials = netrc.netrc()
        login, _, password = credentials.authenticators(self._host)
        return login, password


##############################################################################
class ObjectMissing(Exception):
    """
    Thrown when an expected file object is not found.
    """
    # TODO: Could a standard exception be used here?
    #
    pass  # pylint: disable=unnecessary-pass


##############################################################################
class Uploader(ABC):
    """
    Interface and core functionality for uploading sites.
    """
    @abstractmethod
    def upload(self, file_object, filename):
        """
        Upload the contents of the fileObject to the server as filename.
        """
        raise NotImplementedError()

    @abstractmethod
    def _remove_directory(self, filename):
        """
        Removes the directory, if it exists.
        """
        raise NotImplementedError()

    @abstractmethod
    def _rename_directory(self, original_filename, new_filename):
        """
        If the original filename exists, rename it to the new filename.
        """
        raise NotImplementedError()

    @abstractmethod
    def _make_directory(self, filename):
        """
        Creates a directory on the mirror.
        """
        raise NotImplementedError()

    def prepare(self, directory):
        """
        Prepares to upload to directory.
        """
        # TODO: Should this happen in a constructor? Whatever should happen
        #       accessing these austensibly private members from children
        #       should not.
        #
        self._work_directory = directory + '.new'
        self._final_directory = directory
        self._previous_directory = directory + '-previous'

        self._remove_directory(self._work_directory)
        self._make_directory(self._work_directory)

    def ensure_directory(self, filename):
        """
        Creates a directory on the mirror if none exists.
        """
        self._make_directory(os.path.join(self._work_directory, filename))

    def commit(self):
        """
        Replaces any existing file tree with the uploaded one.
        """
        self._remove_directory(self._previous_directory)
        try:  # self._final_directory may not exist
            self._rename_directory(self._final_directory,
                                   self._previous_directory)
        except ObjectMissing:
            pass
        self._rename_directory(self._work_directory, self._final_directory)


##############################################################################
class UploaderFile(Uploader):
    """
    Uploads to a filesystem.
    """
    def __init__(self):
        pass

    def upload(self, file_object, filename):
        absolute_filename = os.path.join(self._work_directory, filename)
        with open(absolute_filename, 'wb') as destination:
            while True:
                block = file_object.read(1024)
                if block == b'':
                    break
                destination.write(block)
        os.chmod(absolute_filename,
                 stat.S_IRUSR | stat.S_IWUSR | stat.S_IRGRP | stat.S_IROTH)

    def _make_directory(self, filename):
        if not os.path.exists(filename):
            os.mkdir(filename)
            os.chmod(filename,
                     stat.S_IRUSR | stat.S_IWUSR | stat.S_IXUSR
                     | stat.S_IRGRP | stat.S_IXGRP
                     | stat.S_IROTH | stat.S_IXOTH)

    def _remove_directory(self, filename):
        if os.path.exists(filename):
            shutil.rmtree(filename)

    def _rename_directory(self, original_filename, new_filename):
        if os.path.exists(original_filename):
            os.rename(original_filename, new_filename)


##############################################################################
class UploaderFTP(Uploader):
    """
    Uploads to an FTP server.
    """
    def __init__(self, host, credentials=None):
        self._client = ftplib.FTP(host, timeout=30)

        if credentials:
            username, password = credentials.getCredentials()
            self._client.login(username, password)
            del password

    def upload(self, file_object, filename):
        try:
            self._client.storbinary('STOR {}/{}'.format(self._work_directory,
                                                        filename), file_object)
        except Exception as ex:
            message = 'Unable to upload file to server: {}'
            raise Exception(message.format(filename))

    def _make_directory(self, filename):
        objects = self._client.nlst(filename)
        if objects == []:
            self._client.mkd(filename)

    def _remove_directory(self, filename):
        objects = self._client.nlst(filename)
        if objects == []:
            # Object may not exist or may be an empty directory
            try:
                self._client.rmd(filename)
            except ftplib.all_errors:
                # TODO: Clearly we should better understand FTP error handling.
                pass
            return

        if objects == [filename]:
            try:
                self._client.delete(filename)
            except Exception as ex:
                message = 'Unable to delete file "{}" on server: {}'
                raise Exception(message.format(filename, ex))
        else:
            for fobject in objects:
                self._remove_directory(fobject)
            try:
                self._client.rmd(filename)
            except Exception as ex:
                message = 'Unable to delete directory "{}" on server: {}'
                raise Exception(message.format(filename, ex))

    def _rename_directory(self, original_filename, new_filename):
        objects = self._client.nlst(original_filename)
        if objects == []:
            raise ObjectMissing('No such file: {}'.format(original_filename))

        try:
            self._client.rename(original_filename, new_filename)
        except Exception as ex:
            message = 'Unable to rename "{}" as "{}": {}'
            raise Exception(message.format(original_filename,
                                           new_filename,
                                           ex))


##############################################################################
class Mirrorer:  # pylint: disable=too-few-public-methods
    """
    Mirror a web site to a remote server.
    """
    def __init__(self, destination,
                 file_transforms=None, tree_transforms=None):
        """
        Constructor taking upload URL and fileTransforms to apply to HTML
        files.
        """
        url = urlparse(destination)

        ET.register_namespace('', 'http://www.w3.org/1999/xhtml')

        if url.scheme == 'file':
            self._uploader = UploaderFile()
            self._uploader.prepare(url.path)
        elif url.scheme == 'ftp':
            # Slice leading slash off path to make it relative.
            self._uploader = UploaderFTP(url.netloc,
                                         NetrcCredentials(url.netloc))
            self._uploader.prepare(url.path[1:])
        else:
            message = 'Only file and FTP destinations are supported'
            raise NotImplementedError(message)

        if file_transforms is None:
            self._file_transform = []
        else:
            self._file_transform = file_transforms
        if tree_transforms is None:
            self._tree_transform = []
        else:
            self._tree_transform = tree_transforms

    def mirror(self, source):
        """
        Mirror from the provided source.
        """
        ET.register_namespace('', 'http://www.w3.org/1999/xhtml')
        for root, dirs, files in os.walk(source):
            for transformation in self._tree_transform:
                transformation.visit(os.path.relpath(root, source),
                                     dirs,
                                     files)

            for directory in dirs:
                absolute_filename = os.path.join(root, directory)
                relative_filename = os.path.relpath(absolute_filename, source)
                self._uploader.ensure_directory(relative_filename)

            for filename in files:
                absolute_filename = os.path.join(root, filename)
                with open(absolute_filename, 'rb') as file_stream:
                    transformed_stream \
                        = self._transform_file(absolute_filename,
                                               file_stream)
                    self._uploader.upload(transformed_stream,
                                          os.path.relpath(absolute_filename,
                                                          source))
                    transformed_stream.close()

        for transformation in self._tree_transform:
            for filename, content in transformation.newFiles().items():
                fake_file = io.StringIO()
                print(content, file=fake_file)
                fake_file.seek(0, os.SEEK_SET)
                transformed_stream = self._transform_file(filename, fake_file)
                self._uploader.upload(transformed_stream, filename)
                transformed_stream.close()
                fake_file.close()

        self._uploader.commit()

    # This map of common typographical entities to their character
    # representation was derived on 2021-07-21 from
    # https://www.w3.org/wiki/Common_HTML_entities_used_for_typography
    #
    # Although "lt" and "gt" are special in XML they still have to be handled
    # by our manual process.
    #
    _TYPO_ENTITIES = {
        'cent': '#162',
        'pound': '#163',
        'sect': '#167',
        'copy': '#169',
        'laquo': '#171', 'raquo': '#187',
        'reg': '#174',
        'deg': '#176',
        'plusmn': '#177',
        'para': '#182',
        'middot': '#183',
        'frac12': '#188',
        'ndash': '#8211',
        'mdash': '#8212',
        'lsquo': '#8216', 'rsquo': '#8217',
        'sbquo': '#8218',
        'ldquo': '#8220', 'rdquo': '#8221',
        'bdquo': '#8222',
        'dagger': '#8224', 'Dagger': '#8225',
        'bull': '#8226',
        'hellip': '#8230',
        'prime': '#8242', 'Prime': '#8243',
        'euro': '#8364',
        'trade': '#8482',
        'asymp': '#8776',
        'ne': '#8800',
        'le': '#8804', 'ge': '#8805',
        'lt': 'lt', 'gt': 'gt'
        }
    # This map of common entities is built up over time as we gain experience.
    #
    # Although "amp" is special in XML it still has to be handled by our manual
    # process.
    #
    # The "quot" translation results in a literal '"' character which is
    # probably not ideal but is not invalid. Ideally this entity would never
    # be used and the <q>...</q> markup would be used instead. This solution is
    # expedient.
    #
    _COMMON_ENTITIES = {
        'amp': 'amp',
        'nbsp': '#160',
        'quot': '#34'
        }

    def _transform_file(self, filename, stream):
        result = None
        _, extension = os.path.splitext(filename)
        if extension in ['.html', '.xhtml', '.htm']:
            # There used to be a relatively elegant way to do this but Python 3
            # did away with it. Now we have to use this slightly suspect hack
            # of normalising entities by hand.
            #
            def entity_subst(match): return f"&{entity_defs[match.group(1)]};"
            entity_defs = dict(self._COMMON_ENTITIES)
            entity_defs.update(self._TYPO_ENTITIES)
            entity_pattern = re.compile(r'&([^#]+?);')
            processed_doc = []
            try:
                for line in stream:
                    new_line = entity_pattern.sub(entity_subst,
                                                  line.decode('utf-8'))
                    processed_doc.append(new_line)
            except KeyError as ex:
                message = f"Unrecognised entity {ex} found in HTML"
                raise Exception(message)

            try:
                element = ET.fromstringlist(processed_doc)
                tree = ET.ElementTree(element)
            except ET.ParseError as ex:
                message = 'Failed to parse {}: {}'
                raise Exception(message.format(filename, ex))
            tree_root = tree.getroot()

            for transformation in self._file_transform:
                transformation.transform(tree_root)

            content = io.BytesIO()
            tree.write(content, encoding="us-ascii", method="html")
            content.seek(0, os.SEEK_SET)

            # Transform string to file-like 'bytes' object
            stream = content.read()
            result = io.BytesIO(stream)
        else:
            result = stream
        return result
