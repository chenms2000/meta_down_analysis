import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.ResultSet;
import java.sql.Statement;
import java.util.ArrayList;
import java.util.List;

public class BridgeDbDerbyExport {
    private static String clean(String value) {
        if (value == null) {
            return "";
        }
        return value.replace('\t', ' ').replace('\n', ' ').replace('\r', ' ').trim();
    }

    private static String quotedCodeList(String csv) {
        String[] parts = csv.split(",");
        List<String> out = new ArrayList<String>();
        for (String part : parts) {
            String code = part.trim();
            if (!code.matches("[A-Za-z][A-Za-z0-9]*")) {
                throw new IllegalArgumentException("Unsafe BridgeDb code: " + code);
            }
            out.add("'" + code + "'");
        }
        return String.join(",", out);
    }

    private static void exportLinks(Connection conn, String codeCsv) throws Exception {
        String query = "select IDLEFT, CODELEFT, IDRIGHT, CODERIGHT, BRIDGE from APP.LINK";
        if (codeCsv != null && !codeCsv.isBlank()) {
            String codes = quotedCodeList(codeCsv);
            query += " where CODELEFT in (" + codes + ") or CODERIGHT in (" + codes + ")";
        }
        Statement stmt = conn.createStatement();
        stmt.setFetchSize(10000);
        ResultSet rs = stmt.executeQuery(query);
        System.out.println("IDLEFT\tCODELEFT\tIDRIGHT\tCODERIGHT\tBRIDGE");
        while (rs.next()) {
            System.out.println(
                clean(rs.getString(1)) + "\t" +
                clean(rs.getString(2)) + "\t" +
                clean(rs.getString(3)) + "\t" +
                clean(rs.getString(4)) + "\t" +
                clean(rs.getString(5))
            );
        }
        rs.close();
        stmt.close();
    }

    private static void exportInfo(Connection conn) throws Exception {
        Statement stmt = conn.createStatement();
        ResultSet rs = stmt.executeQuery(
            "select SCHEMAVERSION, BUILDDATE, DATASOURCENAME, DATASOURCEVERSION, DATATYPE, SERIES from APP.INFO"
        );
        System.out.println("SCHEMAVERSION\tBUILDDATE\tDATASOURCENAME\tDATASOURCEVERSION\tDATATYPE\tSERIES");
        while (rs.next()) {
            System.out.println(
                clean(rs.getString(1)) + "\t" +
                clean(rs.getString(2)) + "\t" +
                clean(rs.getString(3)) + "\t" +
                clean(rs.getString(4)) + "\t" +
                clean(rs.getString(5)) + "\t" +
                clean(rs.getString(6))
            );
        }
        rs.close();
        stmt.close();
    }

    public static void main(String[] args) throws Exception {
        if (args.length < 2) {
            System.err.println("Usage: BridgeDbDerbyExport <databasePath> <links|info> [codeCsv]");
            System.exit(2);
        }
        Class.forName("org.apache.derby.jdbc.EmbeddedDriver");
        String dbPath = args[0].replace('\\', '/');
        String mode = args[1];
        Connection conn = DriverManager.getConnection("jdbc:derby:" + dbPath);
        conn.setReadOnly(true);
        if ("links".equals(mode)) {
            exportLinks(conn, args.length >= 3 ? args[2] : "");
        } else if ("info".equals(mode)) {
            exportInfo(conn);
        } else {
            throw new IllegalArgumentException("Unknown mode: " + mode);
        }
        conn.close();
        try {
            DriverManager.getConnection("jdbc:derby:" + dbPath + ";shutdown=true");
        } catch (Exception expected) {
            // Derby signals successful shutdown by throwing an exception.
        }
    }
}
