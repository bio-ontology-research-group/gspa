import io
from pathlib import Path
import sys
import tarfile
import tempfile
import unittest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from acquire_gapsmith_archaea import extract


class AcquiredArchiveContract(unittest.TestCase):
    def test_regular_reference_content_is_preserved(self):
        with tempfile.TemporaryDirectory(prefix='gspa-ref-tar-') as directory:
            root=Path(directory);archive=root/'reference.tar.gz'
            with tarfile.open(archive,'w:gz') as tar:
                member=tarfile.TarInfo('rev/1.1.1.1.fasta');content=b'>reference\nMAAA\n'
                member.size=len(content);tar.addfile(member,io.BytesIO(content))
            extract(archive,root/'out')
            self.assertEqual((root/'out/rev/1.1.1.1.fasta').read_bytes(),content)

    def test_traversal_and_link_members_are_rejected_before_extraction(self):
        for kind in ['traversal','link']:
            with self.subTest(kind=kind),tempfile.TemporaryDirectory(prefix='gspa-ref-tar-') as directory:
                root=Path(directory);archive=root/'reference.tar.gz'
                with tarfile.open(archive,'w:gz') as tar:
                    member=tarfile.TarInfo('../outside' if kind=='traversal' else 'link')
                    if kind=='link': member.type=tarfile.SYMTYPE;member.linkname='../outside'
                    tar.addfile(member)
                with self.assertRaises(ValueError): extract(archive,root/'out')
                self.assertFalse((root/'outside').exists())


if __name__=='__main__': unittest.main()
