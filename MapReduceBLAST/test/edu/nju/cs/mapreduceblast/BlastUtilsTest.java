package edu.nju.cs.mapreduceblast;

import static org.junit.Assert.*;

import java.io.IOException;

import org.junit.Test;

public class BlastUtilsTest {

	//@Test
	public void testConvert2SequenceFile() {
		fail("Not yet implemented");
	}

	@Test
	public void testEncode_bp() throws Exception {
		assertTrue((byte)0x00 == BlastUtils.encode_bp('A', true));
		assertTrue((byte)0x01 == BlastUtils.encode_bp('T', false));
		assertTrue((byte)0x10 == BlastUtils.encode_bp('T', true));
		assertTrue((byte)0x20 == BlastUtils.encode_bp('C', true));
		assertTrue((byte)0x30 == BlastUtils.encode_bp('G', true));
		assertTrue((byte)0x40 == BlastUtils.encode_bp('N', true));
	}
	
	@Test (expected=Exception.class)
	public void testEncode_bpInvalidInput() throws Exception{
		BlastUtils.encode_bp('x', true);
	}
	
	@Test
	public void testInt2bytes(){
		byte[] bytes = {(byte) 0xAA,(byte) 0xBB, (byte)0xCC, (byte)0xDD};
		// big-endian
		assertArrayEquals(bytes, BlastUtils.int2bytes(0xAABBCCDD));
	}
	
	@Test
	public void testBytes2int(){
		byte[] bytes = {(byte)0xAA, (byte)0xBB,(byte)0xCC, (byte)0xDD};
//		System.out.printf("types->int: %08X\n", BlastUtils.bytes2int(bytes, 0));
		assertEquals(0xAABBCCDD, BlastUtils.bytes2int(bytes, 0));
	}
	
	@Test
	public void testGetSeqID() throws Exception{
		String seqHeader = ">gi|224589800|ref|NC_000001.10| Homo sapiens " +
				"chromosome 1, GRCh37 primary reference assembly";
		assertEquals("224589800", BlastUtils.getSeqID(seqHeader));
	}
	@Test
	public void testReadBytesFromLoalFile() throws IOException{
		byte[] bytes = BlastUtils.readBytesFromLoalFile("_query.byte");
		assertTrue(bytes.length >0);
	}
}
